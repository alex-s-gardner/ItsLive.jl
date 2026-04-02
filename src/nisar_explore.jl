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