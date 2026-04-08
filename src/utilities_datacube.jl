const datacube_catalog_path = "https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json"

"""
    catalog([catalog_geojson::String])

Return a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes.
Loading is done via GeoDataFrames.jl; when `catalog_geojson` is a URL, it is
fetched with HTTP and then read from a temporary file.

# Examples

```julia
julia> catalog()
```

```julia
julia> catalog("https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json")
```

```julia
julia> catalog("path/to/catalog.json")
```

# Arguments

- `catalog_geojson::String`: URL or local path to the GeoJSON catalog of ITS_LIVE datacubes.

# Author

Alex S. Gardner, JPL, Caltech.
"""
function datacube_catalog(catalog_geojson::String=datacube_catalog_path)
    if startswith(catalog_geojson, "http")
        catalogdf = GDF.read(GDF.ArchGDALDriver(), String(HTTP.get(catalog_geojson).body))
    else
        catalogdf = GDF.read(catalog_geojson)
    end
    return catalogdf
end


function datacube_path(lat, lon; catalog=datacube_catalog_path)
    # Validation
    -90 ≤ lat ≤ 90 || error("lat = $lat, not in range [-90, 90]")
    -180 ≤ lon ≤ 180 || error("lon = $lon, not in range [-180, 180]")

    # Data loading
    pt = GI.Point(lon, lat)
    catalog0 = ItsLive.datacube_catalog(catalog)
    index_datacube = findfirst(GO.intersects.(Ref(pt), catalog0.geometry))

    if isnothing(index_datacube)
        @warn "No intersecting datacube found for point [lat = $lat, lon = $lon]"
        return missing
    end

    path = catalog0[index_datacube, :zarr_url]

    return path
end

function datacube_load(path)
    ds = YAXArrays.open_dataset(Zarr.zopen(path))
    rs = RasterStack((; (k => ds[k] for k in keys(ds.cubes))...); lazy=true)

    return rs
end


function latlon_to_xy(latitude, longitude, rs::AbstractDimStack)
    trans = Proj.Transformation("EPSG:4326", "EPSG:$(DimensionalData.metadata(rs[:mapping])["spatial_epsg"])")
    x, y  = trans(latitude, longitude)
    return x, y
end


# Helper: treat both missing and NaN as invalid
_isvalid(::Missing) = false
_isvalid(x::AbstractFloat) = !isnan(x)
_isvalid(::Any) = true

# Write disaggregation output into a fixed-length destination.
# `disaggregate` ends at max(t2) so its result may be shorter or longer
# than the OutDims grid anchored at output_end. Clip or NaN-pad as needed.
function _write_disagg!(xout::AbstractVector, data::AbstractVector)
    n, m = length(xout), length(data)
    fill!(xout, NaN32)
    xout[1:min(n, m)] .= Float32.(data[1:min(n, m)])
end


# ─── Statistical helpers for dtbias_filter ───────────────────────────────────

# Normal CDF — rational approximation (Abramowitz & Stegun 7.1.26, |ε| < 1.5e-7).
# Avoids any dependency on SpecialFunctions.
function _normcdf(z::Float64)
    t = 1.0 / (1.0 + 0.3275911 * abs(z / sqrt(2.0)))
    p = t * (0.254829592 +
         t * (-0.284496736 +
         t * (1.421413741 +
         t * (-1.453152027 +
         t *  1.061405429))))
    erfc_approx = p * exp(-(z / sqrt(2.0))^2)
    erfc_val = z >= 0.0 ? erfc_approx : 2.0 - erfc_approx
    return 1.0 - erfc_val / 2.0
end

# Tied rank: average-of-ties, O(n log n)
function _tiedrank(x::AbstractVector)
    n  = length(x)
    sp = sortperm(x)
    r  = Vector{Float64}(undef, n)
    i  = 1
    while i <= n
        j = i
        while j < n && x[sp[j+1]] == x[sp[i]]
            j += 1
        end
        avg = (i + j) / 2.0
        for k in i:j
            r[sp[k]] = avg
        end
        i = j + 1
    end
    return r
end

# Spearman one-sided p-value for H₁: r_s < 0 (monotone decrease)
function _spearman_pval(dt::AbstractVector, vp::AbstractVector)
    n = length(dt)
    n < 4 && return 1.0
    rd = _tiedrank(dt)
    rv = _tiedrank(vp)
    rd_mean = sum(rd) / n
    rv_mean = sum(rv) / n
    cov_rv  = sum((rd .- rd_mean) .* (rv .- rv_mean))
    var_d   = sum((rd .- rd_mean).^2)
    var_v   = sum((rv .- rv_mean).^2)
    (var_d == 0.0 || var_v == 0.0) && return 1.0
    r_s   = cov_rv / sqrt(var_d * var_v)
    denom = max(1.0 - r_s^2, 1e-12)
    t     = r_s * sqrt((n - 2) / denom)
    # p-value for H₁: r_s < 0  →  P(Z ≤ t) under normal approximation
    return _normcdf(t)
end

# Mann-Whitney U one-sided p-value for H₁: x stochastically greater than y
function _mannwhitney_pval(x::AbstractVector, y::AbstractVector)
    n1, n2 = length(x), length(y)
    (n1 < 1 || n2 < 1) && return 1.0
    combined = vcat(x, y)
    ranks    = _tiedrank(combined)
    R1 = sum(ranks[1:n1])
    U  = R1 - n1 * (n1 + 1) / 2.0
    mu    = n1 * n2 / 2.0
    sigma = sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    sigma == 0.0 && return 1.0
    z = (U - mu) / sigma
    # p-value for H₁: x > y  →  P(Z ≥ z) = 1 - normcdf(z)
    return 1.0 - _normcdf(z)
end

# Two-sided Mann-Whitney U p-value (bias can go either direction)
function _mannwhitney_pval_twosided(x::AbstractVector, y::AbstractVector)
    p = _mannwhitney_pval(x, y)
    return min(1.0, 2 * min(p, 1.0 - p))
end

# Fisher's method: combine k independent p-values via chi-squared approximation.
# stat = -2 Σ ln(pᵢ) ~ χ²(2k); normal approximation for the CDF tail.
function _fisher_combine(pvals::AbstractVector{Float64})
    isempty(pvals) && return 1.0
    stat = -2.0 * sum(log.(max.(pvals, 1e-300)))
    k    = length(pvals)
    z    = (stat - 2k) / sqrt(4k)
    return 1.0 - _normcdf(z)
end


"""
    sensor_bias_filter(vx, vy, t1, t2, sensor_group_id; ...) → BitVector

Detect and remove observations from sensor groups that produce velocities
significantly slower than the reference sensor. For each sensor group, velocity
is projected onto that group's own mean flow direction, binned by mid-date, and
the co-temporal mean difference vs the reference is tested with a t-like criterion:
a sensor is excluded if `(mean_diff + sescale × SE) < 0`.

Returns a BitVector where `true` = keep observation.

# Arguments
- `vx`, `vy`: 1-D velocity components (m/yr, `Float64`)
- `t1`, `t2`: DateTime vectors (acquisition dates)
- `sensor_group_id`: sensor identifier strings (e.g. from `satellite_img1` field)
- `ref_sensor_group_id`: sensor group id used as reference (default: 1 = Sentinel-2)
- `dtmax_days`: only observations with dt ≤ this value are used for statistics (default: 64.0)
- `bin_width`: temporal bin width as a Period in DateTime space (default: `Month(2)`)
- `mincount`: minimum observations per sensor per bin (default: 3)
- `sescale`: exclusion threshold in standard errors below zero (default: 3.0)

# Author
Alex S. Gardner, JPL, Caltech.
"""
function sensor_bias_filter(
    vx::AbstractVector,
    vy::AbstractVector,
    t1::AbstractVector,
    t2::AbstractVector,
    sensor_group_id::AbstractVector;
    ref_sensor_group_id::Int   = 1,
    dtmax_days::Float64        = 64.0,
    bin_width::Period          = Month(2),
    mincount::Int              = 3,
    sescale::Float64           = 3.0,
)::BitVector

    n    = length(vx)
    keep = trues(n)

    # Early exit: only one (or zero) recognised sensor groups
    unique_ids = unique(sensor_group_id[sensor_group_id .> 0])
    numsg = length(unique_ids)
    length(unique_ids) <= 1 && return keep

    # early exit: reference sensor group not present in data
    !any(unique_ids .== ref_sensor_group_id) && return keep

    # Compute interval (days) and mid-date in DateTime space
    interval_days = Float64.(Dates.value.(t2 .- t1)) ./ 86_400_000.0
    mid_date      = t1 .+ Millisecond.(div.(Dates.value.(t2 .- t1), 2))

    # Select short-dt valid observations for statistics
    valid = (interval_days .<= dtmax_days) .& _isvalid.(vx) .& _isvalid.(vy)
    any(valid) || return keep

    ids_v = sensor_group_id[valid]
    vx_v  = Float64.(vx[valid])
    vy_v  = Float64.(vy[valid])
    md_v  = mid_date[valid]

    # Temporal bin edges anchored to year boundaries in DateTime space
    bin_start = DateTime(year(minimum(md_v)), 1, 1)
    bin_end   = DateTime(year(maximum(md_v)) + 1, 7, 1)
    bin_edges = collect(bin_start:bin_width:bin_end)
    nbins     = length(bin_edges) - 1
    nbins < 1 && return keep

    # Binned mean velocity per sensor group, each projected onto its own mean flow direction
    vbin      = fill(NaN, numsg, nbins)
    vstdbin   = fill(NaN, numsg, nbins)
    vcountbin = zeros(Int, numsg, nbins)

    for i in eachindex(unique_ids)
        sgind = ids_v .== unique_ids[i]
    
        vx0 = Statistics.mean(vx_v[sgind])
        vy0 = Statistics.mean(vy_v[sgind])
        v0  = sqrt(vx0^2 + vy0^2)
        # Skip if v0 is non-finite or below threshold
        (!isfinite(v0) || v0 < 1.0) && continue

        ux = vx0 / v0
        uy = vy0 / v0
        vp    = vx_v[sgind] .* ux .+ vy_v[sgind] .* uy
        md_sg = md_v[sgind]

        for b in 1:nbins
            in_bin = (md_sg .>= bin_edges[b]) .& (md_sg .< bin_edges[b+1])
            count(in_bin) < mincount && continue
            seg = vp[in_bin]
            vcountbin[i, b] = count(in_bin)
            vbin[i, b]      = Statistics.mean(seg)
            vstdbin[i, b]   = length(seg) > 1 ? Statistics.std(seg) : NaN
        end
    end

    # Compare each non-reference group to reference over co-valid bins
    ref_ind = findfirst(unique_ids .== ref_sensor_group_id)
    for i in eachindex(unique_ids)

        unique_ids[i] == ref_sensor_group_id && continue

        covalid = .!isnan.(vbin[ref_ind, :]) .& .!isnan.(vbin[i, :])
        sum(covalid) < 2 && continue

        delta = vbin[i, covalid] .- vbin[ref_ind, covalid]
        m = Statistics.mean(delta)
        s = Statistics.std(delta) / sqrt(sum(covalid) - 1)

        # Sensor is significantly slower than reference → exclude all its observations
        if (m + sescale * s) < 0
            keep[sensor_group_id .== unique_ids[i]] .= false
        end
    end

    return keep
end

"""
    sensor_group_id, sensor_groups = sensor_group(sensor)

return the `sensor` group `sensor_group_id` and the corresponding `sensor_groups`


using Statistics

# Example
```julia
julia> sensor_group_id, sensor_groups = sensor_group(sensor)
```

# Arguments
   - `sensor::Vector{Any}`: sensor list

# Author
Alex S. Gardner, JPL, Caltech.
"""
function sensor_group(sensor)
    sensor_group_id = zeros(Int16, length(sensor))

    for (sg, grp) in sensor_groupings
        for s in grp.sensors
            sensor_group_id[s .== sensor] .= sg
        end
    end

    return sensor_group_id, sensor_groupings
end


# ─── Private helper: MAD-based dt bin filter ─────────────────────────────────
# Returns the maximum dt (days) beyond which the velocity distribution diverges
# from the reference (shortest-dt) bin. Returns Inf when no filtering needed.
function _dtfilter(
    vp::AbstractVector,
    interval_days::AbstractVector{Float64},
    interval_class_edges::Vector{Float64} = [0.0, 16.0, 32.0, 64.0, 128.0, 256.0, 1e10],
    dtbin_mad_thresh::Float64 = 0.5,
)
    min_ref_bin_count = 50

    dt_int = round.(Int, interval_days)
    sp     = sortperm(dt_int)
    dt_s   = dt_int[sp]
    vp_s   = collect(vp[sp])

    bin_ind    = searchsortedlast.(Ref(dt_s), interval_class_edges)
    bin_ind[1] = 1

    nb     = length(interval_class_edges) - 1
    binMed = zeros(Float64, nb)
    binMad = zeros(Float64, nb)
    binCnt = zeros(Int,     nb)

    for i in 1:nb
        lo = max(bin_ind[i], 1)
        hi = bin_ind[i+1]
        lo > hi && continue           # empty bin — no observations in this dt class
        seg = vp_s[lo:hi]
        isempty(seg) && continue
        all(ismissing, seg) && continue
        med        = median(skipmissing(seg))
        binMed[i]  = med
        binMad[i]  = median(skipmissing(abs.(seg .- med)))
        binCnt[i]  = hi - lo + 1
    end

    scale    = dtbin_mad_thresh * 1.4826
    minBound = binMed .- binMad .* scale
    maxBound = binMed .+ binMad .* scale

    ref_ind = findfirst(binCnt .>= min_ref_bin_count)
    isnothing(ref_ind) && (ref_ind = findfirst(maxBound .> 0))
    isnothing(ref_ind) && return Inf

    exclude = (minBound .> maxBound[ref_ind]) .| (maxBound .< minBound[ref_ind])
    return any(exclude) ? Float64(interval_class_edges[findfirst(exclude)]) : Inf
end


"""
    valid_interval, interval_maximum, sensor_groups = interval_bias_filter(vx, vy, interval_millisecond; sensor_group_id)

Identify observations where the velocity distribution shifts with longer image-pair
separation — indicative of "skipping" or "locking" artifacts in feature tracking
glacier velocity estimates.

# Arguments
- `vx`, `vy`: velocity components (m/yr)
- `interval_millisecond`: image-pair separation in **milliseconds**
  (e.g. `Dates.value.(t2 .- t1)`)
- `sensor_group_id`: sensor_group_id identifier; when provided, filtering is applied per group

# Returns
- `valid_interval::BitVector`: `true` = observation passes the filter (keep)
- `interval_maximum`: maximum accepted separation (days) per sensor_group_id group (`Inf` = no limit)
- `sensor_groups`: sensor group metadata from `sensor_group()`

# Author
Alex S. Gardner, JPL, Caltech.
"""
function interval_bias_filter(
    vx::AbstractVector, 
    vy::AbstractVector,
    interval_days::AbstractVector;
    sensor_group_id::Union{AbstractVector,Nothing} = nothing,
    interval_class_edges=[0.0, 16.0, 32.0, 64.0, 128.0, 256.0, Inf],
    min_v0_threshold = 5.0,
    min_count_threshold = 50
)
    n              = length(vx)
    valid_interval = trues(n)

    # Find median flow direction from the shortest-dt observations available
    ind = falses(n)
    for dt0 in interval_class_edges[2:end]
        ind = interval_days .<= dt0
        sum(ind) > min_count_threshold && break
    end

    vx0 = Statistics.median(Float64.(vx[ind]))
    vy0 = Statistics.median(Float64.(vy[ind]))
    v0  = sqrt(vx0^2 + vy0^2)

    # Pixel too slow or non-finite — locking indistinguishable from true signal
    if !isfinite(v0) || v0 < min_v0_threshold
        return valid_interval, Union{Missing,Float64}[Inf], [["none"]]
    end

    ux = vx0 / v0
    uy = vy0 / v0
    vp = Float64.(vx) .* ux .+ Float64.(vy) .* uy

    if isnothing(sensor_group_id)
        interval_maximum_v = _dtfilter(vp, interval_days, interval_class_edges)
        interval_maximum   = Union{Missing,Float64}[interval_maximum_v]
        if interval_maximum_v <= 20_000.0
            valid_interval .&= (interval_days .<= interval_maximum_v) .| ismissing.(vp)
        end
        return valid_interval, interval_maximum, [["none"]]
    else
        unique_ids = unique(sensor_group_id[sensor_group_id.>0])
        interval_maximum = Vector{Union{Missing,Float64}}(missing, length(unique_ids))
        for i in eachindex(unique_ids)
            sgind   = sensor_group_id .== unique_ids[i]
            interval_maximum_v    = _dtfilter(vp[sgind], interval_days[sgind], interval_class_edges)
            interval_maximum[i]  = interval_maximum_v
            if interval_maximum_v <= 20_000.0
                valid_interval .&= (.!sgind) .| (interval_days .<= interval_maximum_v) .| ismissing.(vp)
            end
        end
        return valid_interval, interval_maximum, unique_ids
    end
end





function disaggregate(method, vx, vy, vx_err, vy_err, t1, t2, sensor_group_id;
    output_start=nothing, output_end=nothing,
    output_period=Week(1), loss_norm=TemporalDisaggregations.HuberLoss(1.35),
    sigma_buffer=2, time_buffer=Month(1), verbose=false,
    # IRLS parameters
    irls_max_iter::Int=50,
    irls_tol::Float64=1e-8,
    # Redundancy filtering parameters
    apply_redundancy_filter=false,
    redundancy_interval_bins::Union{Period,AbstractVector{<:Period}}=[Day(0), Day(16), Day(32), Day(64), Day(128), Day(256), Day(1E4)],
    redundancy_temporal_overlap::Float64=0.5,
    redundancy_bin_threshold::Union{Int, Vector{Int}}=20)

    # ── Stage 0: coarse validity screen ──────────────────────────────────────

    # Default output window to the span of valid observations (must happen before computing output_start1/output_end1)
    isnothing(output_start) && (output_start = minimum(t1))
    isnothing(output_end) && (output_end = maximum(t2))

    # Ensure type consistency: if input times are DateTime, convert output times to DateTime
    # This prevents type mismatch in TemporalDisaggregations._date_grid (especially for GP method)
    if eltype(t1) <: DateTime
        output_start = DateTime(output_start)
        output_end = DateTime(output_end)
    end

    output_start1 = output_start - time_buffer
    output_end1 = output_end + time_buffer

    # Fused broadcast for validity checks (reduces allocations)
    valid_obs = @. (abs(vx) < 20000) & (abs(vy) < 200000) &
                   isfinite(vx) & isfinite(vy) &
                   (t2 >= output_start1) & (t1 <= output_end1)
    valid_obs = coalesce.(valid_obs, false)  # Handle missing values
    n_obs1 = sum(valid_obs)

    vx = Float64.(vx[valid_obs])
    vy = Float64.(vy[valid_obs])
    vx_err = Float64.(vx_err[valid_obs])
    vy_err = Float64.(vy_err[valid_obs])
    t1 = t1[valid_obs]
    t2 = t2[valid_obs]
    sensor_group_id = sensor_group_id[valid_obs]

    # Ensure error terms are finite and positive (prevents Inf from 1/err)
    @. vx_err = max(vx_err, 1e-6)  # Replace zeros/negatives with minimum
    @. vy_err = max(vy_err, 1e-6)

    # Filter out non-finite errors (fused broadcast)
    error_valid = @. isfinite(vx_err) & isfinite(vy_err)
    if !all(error_valid)
        vx = vx[error_valid]
        vy = vy[error_valid]
        vx_err = vx_err[error_valid]
        vy_err = vy_err[error_valid]
        t1 = t1[error_valid]
        t2 = t2[error_valid]
        sensor_group_id = sensor_group_id[error_valid]
    end

    # Convert DateTime difference (milliseconds) to fractional days (fused)
    interval_days = @. Float64(Dates.value(t2 - t1)) / 86_400_000.0
    decyear_out = yeardecimal.(output_start1:output_period:output_end1)


    # ── Stage 1: sensor-bias  ────────────────────────────────────────────────
    keep = ItsLive.sensor_bias_filter(vx, vy, t1, t2, sensor_group_id)
    n_obs2 = sum(keep)

    if n_obs2 <= length(decyear_out)
        verbose && @warn "Insufficient observations after sensor filtering: $n_obs2 obs for $(length(decyear_out)) output times"
        return (nothing, nothing, nothing)
    end

    # ── Stage 2: interval-bias filters ───────────────────────────────────────
    keep[keep], _, _ = interval_bias_filter(
        vx[keep], vy[keep], interval_days[keep];
        sensor_group_id=sensor_group_id[keep])
    n_obs3 = sum(keep)

    # ── Validation: Check for finite values before fitting ──────────────────
    # Ensure all values are finite (no NaN/Inf)
    if any(.!isfinite.(vx[keep])) || any(.!isfinite.(vy[keep]))
        return (nothing, nothing, nothing)
    end

    # Ensure error terms are positive (needed for weights = 1/err)
    if any(vx_err[keep] .<= 0) || any(vy_err[keep] .<= 0)
        return (nothing, nothing, nothing)
    end

    # Ensure we have enough observations after filtering
    sum(keep) < 10 && return (nothing, nothing, nothing)

    # ── Stage 3: first fit (build velocity-uncertainty corridor) ─────────────
    # Represent each observation as a horizontal line segment in (time, velocity)
    # space; used later to test geometric intersection with the fit corridor.

    yr_t1 = yeardecimal.(t1)
    yr_t2 = yeardecimal.(t2)

    # Pre-allocate geometry arrays (reduces allocations)
    n_obs = length(yr_t1)
    obs_lines_vx = Vector{GI.Line}(undef, n_obs)
    obs_lines_vy = Vector{GI.Line}(undef, n_obs)
    @inbounds for i in 1:n_obs
        obs_lines_vx[i] = GI.Line([GI.Point(yr_t1[i], vx[i]), GI.Point(yr_t2[i], vx[i])])
        obs_lines_vy[i] = GI.Line([GI.Point(yr_t1[i], vy[i]), GI.Point(yr_t2[i], vy[i])])
    end


    #TODO: impliment observation reduction by excluding redundent interval_values that have similar aquisition times and intervals.
    
    # ── Stage 0.5: redundancy reduction ──────────────────────────────────────
    if apply_redundancy_filter
        # filter is lightning fast

        # Apply interval-stratified redundancy filter (from TemporalDisaggregations.jl)
        keep_vx = copy(keep)
        keep_vx[keep_vx] = TemporalDisaggregations.redundancy_filter(
            vx_err[keep], t1[keep], t2[keep];
            interval_bins=redundancy_interval_bins,
            temporal_overlap=redundancy_temporal_overlap,
            bin_count_threshold=redundancy_bin_threshold)

        keep_vy = copy(keep)
        keep_vy[keep] = TemporalDisaggregations.redundancy_filter(
            vy_err[keep], t1[keep], t2[keep];
            interval_bins=redundancy_interval_bins,
            temporal_overlap=redundancy_temporal_overlap,
            bin_count_threshold=redundancy_bin_threshold)

        if verbose
            pct_kept = round(100 * sum(keep_vx) / n_obs3, digits=1)
            @info "Redundancy filter vx: kept $(sum(keep_vx)) / $n_obs3 observations ($pct_kept% retention)"

            pct_kept = round(100 * sum(keep_vy) / n_obs3, digits=1)
            @info "Redundancy filter vy: kept $(sum(keep_vy)) / $n_obs3 observations ($pct_kept% retention)"

        end
    else
        keep_vx = keep
        keep_vy = keep
    end

    # Fit over extended window (±time_buffer) to reduce edge effects in the corridor.
    # Weights = 1/err give heteroscedastic noise so high-error obs have less influence.
    vx_fit1 = TemporalDisaggregations.disaggregate(
        method, vx[keep_vx], t1[keep_vx], t2[keep_vx];
        output_start=output_start1, output_end=output_end1,
        output_period, loss_norm, irls_max_iter, irls_tol)

    vx_std = std(TemporalDisaggregations.interval_average(vx_fit1, t1[keep], t2[keep]) .- vx[keep])
    vx_line = GI.LineString(GI.Point.(decyear_out, vx_fit1.signal.data))
    vx_line_buff = ItsLive.buffer_unidirectional(vx_line, sigma_buffer * vx_std, dims=2)
    vx_intersects = GO.intersects.(Ref(vx_line_buff), obs_lines_vx)

    vy_fit1 = TemporalDisaggregations.disaggregate(
        method, vy[keep_vy], t1[keep_vy], t2[keep_vy];
        output_start=output_start1, output_end=output_end1,
        output_period, loss_norm, irls_max_iter, irls_tol)

    #valid_obs[valid_obs] = keep
    #return (vx_fit1, vy_fit1, valid_obs)

    vy_std = std(TemporalDisaggregations.interval_average(vy_fit1, t1[keep], t2[keep]) .- vy[keep])
    vy_line = GI.LineString(GI.Point.(decyear_out, vy_fit1.signal.data))
    vy_line_buff = ItsLive.buffer_unidirectional(vy_line, sigma_buffer * vy_std, dims=2)
    vy_intersects = GO.intersects.(Ref(vy_line_buff), obs_lines_vy)

    # ── Stage 3: corridor re-inclusion ───────────────────────────────────────
    # Re-include any observation (even Stage-1 rejects) that intersects both corridors.
    keep = (vx_intersects .& vy_intersects)

    # ── Stage 0.5: redundancy reduction ──────────────────────────────────────
    if apply_redundancy_filter
        # filter is lightning fast

        # Apply interval-stratified redundancy filter (from TemporalDisaggregations.jl)
        keep_vx = copy(keep)
        keep_vx[keep_vx] = TemporalDisaggregations.redundancy_filter(
            vx_err[keep], t1[keep], t2[keep];
            interval_bins=redundancy_interval_bins,
            temporal_overlap=redundancy_temporal_overlap,
            bin_count_threshold=redundancy_bin_threshold)

        keep_vy = copy(keep)
        keep_vy[keep] = TemporalDisaggregations.redundancy_filter(
            vy_err[keep], t1[keep], t2[keep];
            interval_bins=redundancy_interval_bins,
            temporal_overlap=redundancy_temporal_overlap,
            bin_count_threshold=redundancy_bin_threshold)

        if verbose
            pct_kept = round(100 * sum(keep_vx) / n_obs3, digits=1)
            @info "Redundancy filter vx: kept $(sum(keep_vx)) / $n_obs3 observations ($pct_kept% retention)"

            pct_kept = round(100 * sum(keep_vy) / n_obs3, digits=1)
            @info "Redundancy filter vy: kept $(sum(keep_vy)) / $n_obs3 observations ($pct_kept% retention)"

        end
    else
        keep_vx = keep
        keep_vy = keep
    end

    # ── Stage 4: final fit on expanded observation set ───────────────────────
    # Pre-allocate weight arrays (reduces allocations)
    weights_vx = similar(vx_err, sum(keep_vx))
    @inbounds @. weights_vx = 1 / vx_err[keep_vx]^2

    weights_vy = similar(vy_err, sum(keep_vy))
    @inbounds @. weights_vy = 1 / vy_err[keep_vy]^2

    vx_fit = TemporalDisaggregations.disaggregate(
        method, vx[keep_vx], t1[keep_vx], t2[keep_vx];
        output_start, output_end, output_period, loss_norm, irls_max_iter, irls_tol, weights=weights_vx)

    vy_fit = TemporalDisaggregations.disaggregate(
        method, vy[keep_vy], t1[keep_vy], t2[keep_vy];
        output_start, output_end, output_period, loss_norm, irls_max_iter, irls_tol, weights=weights_vy)

    # Map keep indices back to the original (pre-Stage-0) index space
    valid_obs[valid_obs] = keep

    if verbose
        pct_sensor = round(Int, (n_obs1 - n_obs2) / n_obs1 * 100)
        pct_interval = round(Int, (n_obs2 - n_obs3) / n_obs1 * 100)
        pct_sigma = -round(Int, (n_obs3 - sum(keep)) / n_obs1 * 100)
        println(rpad("Filter", 20), rpad("% removed", 10))
        println("-"^30)
        println(rpad("sensor_bias", 20), pct_sensor)
        println(rpad("interval_bias", 20), pct_interval)
        println(rpad("sigma_envelope", 20), pct_sigma)
        println()

    end
    return (vx_fit, vy_fit, valid_obs)
end
