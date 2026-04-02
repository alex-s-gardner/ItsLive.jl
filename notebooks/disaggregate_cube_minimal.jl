
# ── 1. Dependencies ──────────────────────────────────────────────────────────
begin
    using ItsLive
    using TemporalDisaggregations
    using YAXArrays
    using Zarr
    using DimensionalData
    using Dates
    using LinearAlgebra
    using CairoMakie
end

# Zarr.MaxLengthStrings.length() returns `nothing` for all-zero (empty) fill-value strings,
# which propagates into ncodeunits and breaks JSON serialization via `1:nothing`.
# Fix: return 0 instead of nothing for empty MaxLengthStrings.
Base.ncodeunits(s::Zarr.MaxLengthStrings.MaxLengthString) = something(findlast(!iszero, s.data), 0)

# ── 2. Locate datacube ───────────────────────────────────────────────────────
begin
    lat, lon = 60.0626, -139.3193   # Hubbard Glacier, Alaska
    zarr_path = ItsLive.datacube_path(lat, lon)
    println("Zarr path: ", zarr_path)
end

ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path))
ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path), skip_keys=setdiff(keys(ds.cubes), (:vx, :vy,:acquisition_date_img1, :acquisition_date_img2, :satellite_img1)))

local_dir = "/Users/gardnera/data/its-live-data/datacubes/v02"
zarr_path_local = joinpath(local_dir, last(splitpath(zarr_path)))
#savedataset(ds; path=zarr_path_local, driver=:zarr, overwrite=true)


# ── 3. Configure disaggregation ──────────────────────────────────────────────
begin
    method = Spline(smoothness=1e-3, tension=1)
    output_dir = joinpath(tempdir(), "itslive_disagg")
    output_start = Date(2014, 1, 1)
    output_end = Date(2024, 1, 1)
    output_period = Week(1)
end

vars = [:vx, :vy]

# Capture date metadata before entering mapCube (closure capture)
t1 = collect(ItsLive.to_datetime(ds["acquisition_date_img1"]))
t2 = collect(ItsLive.to_datetime(ds["acquisition_date_img2"]))
sensor = collect(ds["satellite_img1"])

# Build output time axis
output_times = collect(output_start:output_period:output_end)

mkpath(output_dir)


output_start = Date(2014, 1, 1)
output_end = Date(2025, 1, 1)
output_period = Week(1)


function f2(vx_col, vy_col)
    method = Spline(smoothness=1e-3, tension=10)
    output_start = Date(2014, 1, 1)
    output_end = Date(2025, 1, 1)
    output_period = Week(1)
    loss_norm = :L1
    max_abs_velocity = 20000.0
    min_valid_observations = 10

    f = Base.Fix{3}(Base.Fix{4}(Base.Fix{5}(ItsLive.disaggregate, sensor), t2), t1)

    result = try
        f(vx_col, vy_col; method, output_start, output_end, output_period, loss_norm, max_abs_velocity, min_valid_observations)
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        nothing
    end

    if isnothing(result)
        n = length(output_start:output_period:output_end)
        return fill(NaN,n)
    else
        return result.signal.data
    end
end

vx0 = ds.vx[700:end, 800:end, :];
vy0 = ds.vy[700:end, 800:end, :];

Ti_out = Ti(output_start:output_period:output_end)

out0 = XOutput(Ti_out; outtype=Float64)

r = xmap(f2, vx0 ⊘ :mid_date, vy0 ⊘ :mid_date, output=out0, inplace=false)

gen_cube = compute_to_zarr(Dataset(layer=r), "my_gen_cube.zarr", overwrite=true, max_cache=1e9)

