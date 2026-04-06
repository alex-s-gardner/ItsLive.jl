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
    using MultivariateStats

    # For faster BLAS on macOS: 
    using AppleAccelerate
    using Profile
end


# ── 2. Locate datacube ───────────────────────────────────────────────────────
begin
    zarr_path = "https://its-live-data.s3.amazonaws.com/test-space/datacubes/nisar/datacubes-04.01.2026/N60W130/ITS_LIVE_vel_EPSG3413_G0120_X-3250000_Y250000.zarr"

    ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path))
    ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path), skip_keys=setdiff(keys(ds.cubes), (:vx, :vy, :vx_error, :vy_error, :acquisition_date_img1, :acquisition_date_img2, :satellite_img1)))
    #ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path), skip_keys=setdiff(keys(ds.cubes), (:vx,)))


    local_dir = "/Users/gardnera/data/its-live-data/datacubes/v02"
    #local_dir = "/mnt/devon-r3/data/test"

    zarr_path_local = joinpath(local_dir, last(splitpath(zarr_path)))
    #savedataset(ds; path=zarr_path_local, overwrite=true)
    # Bash(aws s3 sync --no-sign-request s3://its-live-data/test-space/datacubes/nisar/datacubes-04.01.2026/N60W130/ITS_LIVE_vel_EPSG3413_G0120_X-3250000_Y250000.zarr /Users/gardnera/data/its-live-data/datacubes/v02/)
    
    # Discover the datacube zarr URL from the ITS_LIVE catalog
    rs     = ItsLive.datacube_load(zarr_path)
end;


#create plots to visualize the data, and apply filters to extract meaningful information from the datacube.

granule_url = collect(rs[:granule_url])
rslc = occursin.("RSLC", granule_url);

# Create histogram plot comparing RSLC and GSLC error distributions
bins = 0:10:250
begin
    vx_error = collect(rs[:vx_error])
    vy_error = collect(rs[:vy_error])
    gslc = .!rslc

    # Filter out missing values
    vx_rslc = vx_error[rslc .& .!ismissing.(vx_error)]
    vx_gslc = vx_error[gslc .& .!ismissing.(vx_error)]
    vy_rslc = vy_error[rslc .& .!ismissing.(vy_error)]
    vy_gslc = vy_error[gslc .& .!ismissing.(vy_error)]

    fig = Figure(size=(1200, 500))

    # Left panel: vx_error
    ax1 = Axis(fig[1, 1],
        xlabel="vx error (m/yr)",
        ylabel="Count",
        title="vx std over stable surface")
    hist!(ax1, vx_rslc; color=(:blue, 0.5), label="RSLC [n = $(length(vx_rslc))]", bins)
    hist!(ax1, vx_gslc; color=(:red, 0.5), label="GSLC [n = $(length(vx_gslc))]", bins)
    axislegend(ax1, position=:rt)

    # Right panel: vy_error
    ax2 = Axis(fig[1, 2],
        xlabel="vy error (m/yr)",
        ylabel="Count",
        title="vy std over stable surface")
    hist!(ax2, vy_rslc; color=(:blue, 0.5), label="RSLC [n = $(length(vy_rslc))]", bins)
    hist!(ax2, vy_gslc; color=(:red, 0.5), label="GSLC [n = $(length(vy_gslc))]", bins)
    axislegend(ax2, position=:rt)

    fig
end

vx = rs[:vx]
vy = rs[:vy]

x = dims(rs[:vx], :x)
y = dims(rs[:vx], :y)

# ── 3. Disaggregation configuration ─────────────────────────────────────────
begin
    output_period = Month(1)
    loss_norm = :L1
    sigma_buffer = 2
    time_buffer = Month(1)
    observation_error_minimum = 20 # m/yr
    verbose = false

    method = Spline(smoothness=1e-3, tension=1)

end


begin
    t1 = collect(ItsLive.to_datetime(rs["acquisition_date_img1"]))
    t2 = collect(ItsLive.to_datetime(rs["acquisition_date_img2"]))

    vx_err = collect(ds[:vx_error])
    vy_err = collect(ds[:vy_error])

    sensor = String.(collect(rs["satellite_img1"]))
    sensor_group_id, sensor_groups = ItsLive.sensor_group(sensor)

    # mission_all = collect(rs["mission_img1"]).. error reading mission

    # set minimum error
    vx_err[vx_err.<observation_error_minimum] .= observation_error_minimum
    vy_err[vy_err.<observation_error_minimum] .= observation_error_minimum

    # sort data by ts
    sorted = true
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
    output_start = Date(minimum(skipmissing(t2))) - time_buffer
    output_end = Date(maximum(skipmissing(t2))) - time_buffer

    # Coarse validity screen
    time_output = output_start:output_period:output_end

    m = length(time_output)
    n, p, _ = size(rs[:vx])
end;


# degine output array
begin
    x = X(dims(rs[:vx], :x).val.data)
    y = Y(dims(rs[:vx], :y).val.data)
    ti = Ti(time_output)
    vx_out = zeros(ti, x, y)
    vy_out = zeros(ti, x, y)

    vx_out = zeros(ti, x, y)
    vy_out = zeros(ti, x, y)

    vx_out_data = vx_out.data
    vy_out_data = vy_out.data

end;

@showprogress desc = "Fitting splines to cube..." for i in eachindex(x)[162:end]

    # --------- Pre-load the full (y × time) slice into memory before the threaded loop.
    # X(i) is positional indexing (i-th element along x); avoids per-pixel Zarr reads
    # inside @threads and ensures each Zarr chunk is read at most once per x step.

    # WARNING: this may cause out-of-memory errors if the datacube is too large or the 
    # chunking is inefficient. In that case, consider reading smaller y-slices or individual 
    # pixels inside the loop, at the cost of more I/O.

    vx_col = permutedims(collect(rs[:vx][x=i]))::Array{Union{Missing,Int16},2}  # (nt × ny): each column is one pixel's time series
    vy_col = all(ismissing.(vx_col)) ? continue : permutedims(collect(rs[:vy][x=i]))::Array{Union{Missing,Int16},2}

    @time Threads.@threads :greedy for j in eachindex(y)  # loop over y pixels; adjust range as needed for testing

        all(ismissing.(vy_col[:, j])) && continue  # skip if all vy values are missing (assumes vx is also missing, but check just in case)

        try
            vx_fit, vy_fit, valid_obs = ItsLive.disaggregate(
                method, vx_col[:, j], vy_col[:, j], vx_err, vy_err, t1, t2, sensor_group_id;
                output_start, output_end, output_period, loss_norm, sigma_buffer, time_buffer, verbose)

            if !isnothing(vx_fit)
                vx_out_data[:, i, j] = vx_fit.signal.data
                vy_out_data[:, i, j] = vy_fit.signal.data
            end
        catch e
            println("Error at pixel (x=$i; y=$j) ------------------------------------------")
            rethrow(e)
        end
    end
end




include("nisar_list_files.jl")
foo = list_nisar_nc_files();
foo = s3_to_https.(foo);


# Location info
begin
    latlon = (60.0828, -140.4767)
    pt_name = "Seward Glacier, Alaska"
    detail_interval = (Date(2020, 1, 1), Date(2023, 1, 1))

    latlon = (61.4471, -140.6709)
    pt_name = "Glacier?, Alaska"
    detail_interval = (Date(2016, 1, 1), Date(2019, 1, 1))

    latlon = (60.0510, -139.3518)
    pt_name = "Hubbard Glacier, Alaska"
    detail_interval = (Date(2016, 1, 1), Date(2019, 1, 1))
end

# load datacube
begin
    path = ItsLive.datacube_path(latlon...)
    rs = ItsLive.datacube_load(path)

    # find x y and projected coordinates
    x, y = ItsLive.latlon_to_xy(latlon..., rs)

    aggregate_values = rs[:v][x=Near(x), y=Near(y)]
    t1 = ItsLive.to_datetime(rs["acquisition_date_img1"])
    t2 = ItsLive.to_datetime(rs["acquisition_date_img2"])

    valid_index = .!ismissing.(aggregate_values)
    aggregate_values = aggregate_values[valid_index]
    t1 = t1[valid_index]
    t2 = t2[valid_index]
end


