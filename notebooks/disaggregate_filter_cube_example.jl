# # disaggregate_filter_example.jl
#
# Demonstrates ItsLive single-pixel filtering and disaggregation with:
#   • sensor-bias filter  (ItsLive.sensor_bias_filter)
#   • dt-bias filter      (ItsLive.dtbias_filter)
#   • two-stage re-inclusion of valid long-dt observations
#   • single-pixel plot   (raw data + filter mask + disaggregated signals)
#
# Run cell-by-cell in VS Code (Julia extension ≥ 1.79) or top-to-bottom as a
# plain script.  A working internet connection is required to stream the
# ITS_LIVE S3 datacube.
#
# Author: Alex S. Gardner, JPL / Caltech.

# ── 1. Dependencies ─────────────────────────────────────────────────────────
begin
    using ItsLive
    using TemporalDisaggregations
    using YAXArrays
    using Zarr
    using DimensionalData
    using Dates
    using Statistics
    using CairoMakie
    using Colors
    import GeoInterface as GI
    import GeometryOps as GO
    using DateFormats: yeardecimal
    using KernelFunctions
    using LinearAlgebra
    using ProgressMeter

    # For faster BLAS on macOS: 
    using AppleAccelerate 
    using Profile
            
    #For faster BLAS on Intel Linux: 
    #using MKL

    # Zarr.MaxLengthStrings.length() returns `nothing` for all-zero (empty) fill-value strings,
    # which propagates into ncodeunits and breaks JSON serialization via `1:nothing`.
    # Fix: return 0 instead of nothing for empty MaxLengthStrings.
    Base.ncodeunits(s::Zarr.MaxLengthStrings.MaxLengthString) = something(findlast(!iszero, s.data), 0)
end

# ── 3. Disaggregation configuration ─────────────────────────────────────────
begin
    output_period = Week(1)
    loss_norm = :L1
    sigma_buffer = 2
    time_buffer = Month(1)
    observation_error_minimum = 20 # m/yr
    verbose = false

    method = Spline(smoothness=1e-3, tension=1)

end


# ── 2. Locate datacube ───────────────────────────────────────────────────────
begin
    lat, lon = 60.0626, -139.3193   # Hubbard Glacier, Alaska
    zarr_path = ItsLive.datacube_path(lat, lon)
end

begin
    ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path))
    ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path), skip_keys=setdiff(keys(ds.cubes), (:vx, :vy, :vx_error, :vy_error, :acquisition_date_img1, :acquisition_date_img2, :satellite_img1)))

    local_dir = "/Users/gardnera/data/its-live-data/datacubes/v02"
    #local_dir = "/mnt/devon-r3/data/test"

    zarr_path_local = joinpath(local_dir, last(splitpath(zarr_path)))
    #savedataset(ds; path=zarr_path_local, driver=:zarr, overwrite=true)

    # Discover the datacube zarr URL from the ITS_LIVE catalog
    rs     = ItsLive.datacube_load(zarr_path_local)
    #rs     = ItsLive.datacube_load(zarr_path)

end;

# ── 4. Extract pixel data & apply filters ───────────────────────────────────

# extract raw pixel values and observation intervals
begin
    t1 = collect(ItsLive.to_datetime(rs["acquisition_date_img1"]))
    t2 = collect(ItsLive.to_datetime(rs["acquisition_date_img2"]))

    vx_err = collect(ds[:vx_error])
    vy_err = collect(ds[:vy_error])

    sensor =  String.(collect(rs["satellite_img1"]))
    sensor_group_id, sensor_groups = ItsLive.sensor_group(sensor)
    
    # mission_all = collect(rs["mission_img1"]).. error reading mission

    # set minimum error
    vx_err[vx_err .< observation_error_minimum] .= observation_error_minimum
    vy_err[vy_err .< observation_error_minimum] .= observation_error_minimum

    # sort data by ts
    sorted = true;
    if !issorted(t1)
        sorted = false
        order = sortperm(t1)
        t1 = t1[order]
        t2 = t2[order]
        vx_err = vx_err[order]
        vy_err = vy_err[order]
        sensor = sensor[order]
        
    end

    # define output start and end
    output_start = Date(2014,1,1) + time_buffer
    output_end = Date(maximum(skipmissing(t2))) - time_buffer

    # Coarse validity screen
    time_output = output_start:output_period:output_end

    m = length(time_output)
    n, p,_ = size(rs[:vx])
end;

# degine output array
begin
    x = X(dims(rs[:vx],:x).val.data);
    y = Y(dims(rs[:vx],:y).val.data);
    ti = Ti(time_output);
    vx_out = zeros(ti, x, y);
    vy_out = zeros(ti, x, y);

    vx_out = zeros(ti, x, y);
    vy_out = zeros(ti, x, y);

    vx_out_data = vx_out.data;
    vy_out_data = vy_out.data;

end;


@showprogress desc = "Fitting splines to cube..." for i in eachindex(x)[162:end]
#i = 258
    # --------- Pre-load the full (y × time) slice into memory before the threaded loop.
    # X(i) is positional indexing (i-th element along x); avoids per-pixel Zarr reads
    # inside @threads and ensures each Zarr chunk is read at most once per x step.
    
    # WARNING: this may cause out-of-memory errors if the datacube is too large or the 
    # chunking is inefficient. In that case, consider reading smaller y-slices or individual 
    # pixels inside the loop, at the cost of more I/O.

    vx_col = permutedims(collect(rs[:vx][x=i]))::Array{Union{Missing, Int16}, 2}  # (nt × ny): each column is one pixel's time series
    vy_col = permutedims(collect(rs[:vy][x=i]))::Array{Union{Missing, Int16}, 2}

    # Optional: filter observations by granule type
    granule_url = collect(rs[:granule_url])
    rslc = occursin.("RSLC", granule_url)
    gslc = .!rslc
    # use_index = gslc  # Uncomment to use only GSLC observations
    use_index = trues(length(t1))  # Use all observations

    # For single pixel testing, uncomment and set j:
    j = 400
    all(ismissing.(vy_col[:, j])) && error("All missing values at pixel j=$j")

    try
        vx_fit, vy_fit, valid_obs = ItsLive.disaggregate(
            method, vx_col[use_index, j], vy_col[use_index, j], vx_err[use_index], vy_err[use_index],
            t1[use_index], t2[use_index], sensor_group_id[use_index];
            output_start, output_end, output_period, loss_norm, sigma_buffer, time_buffer, verbose)

        if !isnothing(vx_fit)
            vx_out_data[:,i,j] = vx_fit.signal.data
            vy_out_data[:,i,j] = vy_fit.signal.data
        end
    catch e
        println("Error at pixel (x=$i; y=$j) ------------------------------------------")
        rethrow(e)
    end
end






# ── 5. Figure — single-pixel: raw data, removed points, disaggregated signalbegin
    fig1 = ItsLive.plot_disaggregate_pixel(
        vx, vy, t1, t2, vx_fit, vy_fit, valid_obs;
        method, pt_name, latlon, detail_interval)

    display(fig1)
    save(joinpath(@__DIR__, "pixel_disagg_comparison.png"), fig1, px_per_unit=2)
    println("Saved: pixel_disagg_comparison.png")
end
end
