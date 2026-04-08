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

    # For faster BLAS on macOS: 
    using AppleAccelerate 
            
    #For faster BLAS on Intel Linux: 
    #using MKL

    # Zarr.MaxLengthStrings.length() returns `nothing` for all-zero (empty) fill-value strings,
    # which propagates into ncodeunits and breaks JSON serialization via `1:nothing`.
    # Fix: return 0 instead of nothing for empty MaxLengthStrings.
    Base.ncodeunits(s::Zarr.MaxLengthStrings.MaxLengthString) = something(findlast(!iszero, s.data), 0)
end


# ── 3. Disaggregation configuration ─────────────────────────────────────────
begin
    output_period = Month(1)
   
    loss_norm      = TemporalDisaggregations.HuberLoss(1.35)
    #loss_norm      = TemporalDisaggregations.L1DistLoss()
    #loss_norm      = TemporalDisaggregations.L2DistLoss()

    sigma_buffer = 2
    time_buffer = Month(1)
    observation_error_minimum = 20 # m/yr
    verbose = false

    method_selection = :Spline

    if method_selection == :GP
           k =
               5.0^2 * PeriodicKernel(r=[0.5]) * with_lengthscale(Matern52Kernel(), 0.3) +
                50.0^2 * with_lengthscale(Matern52Kernel(), 3.0)

           # k = 100.0^2 * with_lengthscale(Matern52Kernel(), 3.0)

            method = GP(
                kernel=k,
                obs_noise=1.0 ^2,
            )

    elseif method_selection == :Spline
        method = Spline(smoothness=0.25, tension=0.5)
    elseif method_selection == :Sinusoid
    method = Sinusoid()
    end
end

for location in 1:5
#location = 5
begin
    if location == 1
        latlon          = (60.2522, -141.1039)
        pt_name         = "Skipping Example Glacier, Alaska"
        detail_interval = (Date(2020, 6, 1), Date(2022, 6, 1))
    elseif location == 2
        latlon = (-75.181452, -105.993140)
        pt_name = "E. Twaites"
        detail_interval = (Date(2023, 6, 1), Date(2025, 1, 1))
    elseif location == 3
        latlon = (69.1043, -49.4651)
        pt_name = "Jakobshavn"
        detail_interval = (Date(2014, 1, 1), Date(2016, 1, 1))
    elseif location == 4
        latlon = (60.0249, -140.1264)
        pt_name = "Alaska, rock"
        detail_interval = (Date(2014, 1, 1), Date(2016, 1, 1))
    elseif location == 5
        latlon = (60.0855, -140.1812)
        pt_name = "Alaska, surge"
        detail_interval = (Date(2020, 6, 1), Date(2023, 6, 1))
    end

    # Discover the datacube zarr URL from the ITS_LIVE catalog
    path   = ItsLive.datacube_path(latlon...)
    rs     = ItsLive.datacube_load(path)
    x0, y0 = ItsLive.latlon_to_xy(latlon..., rs)
end;

# ── 4. Extract pixel data & apply filters ───────────────────────────────────
begin
    # extract raw pixel values and observation intervals
    vx_err = Float64.(collect(rs[:vx_error]))
    vy_err = Float64.(collect(rs[:vy_error]))

    t1 = collect(ItsLive.to_datetime(rs["acquisition_date_img1"]))
    t2 = collect(ItsLive.to_datetime(rs["acquisition_date_img2"]))

    sensor =  String.(collect(rs["satellite_img1"]))
    sensor_group_id, sensor_groups = ItsLive.sensor_group(sensor)
    # mission_all = collect(rs["mission_img1"]).. error reading mission

    # set minimum error
    vx_err[vx_err .< observation_error_minimum] .= observation_error_minimum
    vy_err[vy_err .< observation_error_minimum] .= observation_error_minimum


    vx = collect(rs[:vx][x=Near(x0), y=Near(y0)])::Vector{Union{Missing,Int16}}
    vy = collect(rs[:vy][x=Near(x0), y=Near(y0)])::Vector{Union{Missing,Int16}}

    # define output start and end
    output_start = Date(2014,1,1) + time_buffer
    output_end = Date(maximum(skipmissing(t2))) - time_buffer


    # Coarse validity screen
    @time (vx_fit, vy_fit, valid_obs) = ItsLive.disaggregate(method, vx, vy, vx_err, vy_err, t1, t2, sensor_group_id; output_start, output_end, output_period, loss_norm, sigma_buffer, time_buffer, verbose=false, apply_redundancy_filter=false, irls_max_iter = 10, irls_tol = 1e-5);
end;

# ── 5. Figure — single-pixel: raw data, removed points, disaggregated signals
begin
    fig1 = ItsLive.plot_disaggregate_pixel(
        vx, vy, t1, t2, vx_fit, vy_fit, valid_obs;
        method, pt_name, latlon, detail_interval)

    display(fig1)
end
end
