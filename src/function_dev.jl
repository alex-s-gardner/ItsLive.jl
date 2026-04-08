begin
    using ItsLive
    using GeoDataFrames
    using GeoInterface
    using GeometryOps
    using Rasters
    using YAXArrays
    using Zarr
    using DimensionalData
    using CairoMakie
    using GeometryBasics
    using DateFormats
    import DimensionalData as DD
    using TemporalDisaggregations
    using Proj
    using KernelFunctions
    using Colors
end


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


    latlon = (-75.181452, -105.993140)
    pt_name = "E. Twaites"
    detail_interval = (Date(2025, 1, 1), Date(2026, 1, 1))
end

# load datacube
begin
    path = ItsLive.datacube_path(latlon...)
    rs = ItsLive.datacube_load(path)

    # find x y and projected coordinates
    x, y = ItsLive.latlon_to_xy(latlon..., rs)

    aggregate_values = rs[:v][x=Near(x), y=Near(y)]
    t1  = ItsLive.to_datetime(rs["acquisition_date_img1"])
    t2  = ItsLive.to_datetime(rs["acquisition_date_img2"])

    valid_index  = .!ismissing.(aggregate_values)
    aggregate_values    = aggregate_values[valid_index]
    t1                  = t1[valid_index]
    t2                  = t2[valid_index]
end

#valid_index = abs.(aggregate_values) .> 200
#aggregate_values = aggregate_values[valid_index]
#t1 = t1[valid_index]
#t2 = t2[valid_index]


# disaggregation
begin
    if false
        methods = (
            Spline(
                smoothness=1e-100 ,
                tension=0
                ),
        )

    else
        # specify GP kernel
        #=
        ┌─────────────────────────────────────────────────────┬───────────────┬──────────────────────────────────────────────────────┐
        │                        Term                         │ Amplitude (σ) │                       Meaning                        │
        ├─────────────────────────────────────────────────────┼───────────────┼──────────────────────────────────────────────────────┤
        │ 15.0^2 * PeriodicKernel * Matern52(lengthscale=3yr) │ 15            │ Annual seasonality that slowly evolves over ~3 years │
        ├─────────────────────────────────────────────────────┼───────────────┼──────────────────────────────────────────────────────┤
        │ 5.0^2 * Matern52(lengthscale=2yr)                   │ 5             │ Smooth multi-year trend                              │
        ├─────────────────────────────────────────────────────┼───────────────┼──────────────────────────────────────────────────────┤
        │ 3.0^2 * Matern32(lengthscale=1/12yr)                │ 3             │ Short-term noise / sub-monthly variation             │
        └─────────────────────────────────────────────────────┴───────────────┴──────────────────────────────────────────────────────┘
        =#
        k =
            2.0^2 * PeriodicKernel(r=[0.5]) * with_lengthscale(Matern52Kernel(), 3.0) +
            15.0^2 * with_lengthscale(Matern52Kernel(), 2.0) +
            10.0^2 * with_lengthscale(Matern32Kernel(), 1 / 12)

        # specify disaggregation methods to run
        methods = (
            Spline(
                smoothness=1e-3 ,
                tension = 1
                ),

            GP(
                kernel   = k,
                obs_noise = 25.0,
                n_quad    = 3
                ),

            Sinusoid()
        )
    end

    # specify disaggregation output parameters
    output_start    = Date(2014, 1, 1)
    output_period   = Week(1)
    loss_norm       = :L2

    # dissaggregate 
    results = TemporalDisaggregations.disaggregate.(
        methods,
        Ref(aggregate_values), Ref(t1), Ref(t2); 
        output_period, output_start, loss_norm
        )
end

# plot
begin
    title_suffex = ": $(pt_name) [$(round(latlon[1], digits=2))°, $(round(latlon[2], digits=2))°]"

    show_legend     = false
    ylabel = "velocity (m yr⁻¹)"

    fig = ItsLive.plot_disaggregate(aggregate_values, t1, t2, results;
        title_suffex,
        detail_interval,
        show_legend,
        ylabel,
    )
end

# save
save("$(pt_name)_disaggregation.png", fig, px_per_unit=2)


path = "https://its-live-data.s3.amazonaws.com/test-space/datacubes/nisar/datacubes/N60W050/ITS_LIVE_vel_EPSG3413_G0120_X-250000_Y-2450000.zarr"
ds = YAXArrays.open_dataset(Zarr.zopen(path))
rs = RasterStack((; (k => ds[k] for k in keys(ds.cubes))...); lazy=true)
valid = (rs.vx_error[:] .< 200) .& (rs.vy_error[:] .< 200)


latlon = (60.0828, -140.4767)
pt_name = "Seward Glacier, Alaska"
detail_interval = (Date(2020, 1, 1), Date(2023, 1, 1))
(x,y) = ItsLive.latlon_to_xy(latlon..., rs)

rs = rs[x=Near(x), y=Near(y)]


v = collect(rs[:v][:, :, 1])
index1 = .!ismissing.(v)
index2 = abs.(v[index1]) .< 2000

fig, ax, hm = heatmap(v; colorrange = (0, 100))
fig.attributes.title = "layer 1 velocity [m/yr]"
Colorbar(fig[:, end+1], hm)