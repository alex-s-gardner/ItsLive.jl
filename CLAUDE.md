# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

ItsLive.jl is a Julia package for accessing, filtering, and analyzing NASA ITS_LIVE satellite glacier velocity and elevation change data. It provides cloud-native access to zarr datacubes on S3, statistical filtering of satellite-derived velocity data, temporal disaggregation, and visualization tools.

## Development Commands

```julia
# Start Julia with the project environment
julia --project=.

# Instantiate dependencies
using Pkg; Pkg.instantiate()

# Load the package in development
using Pkg; Pkg.develop(path=".")
using ItsLive

# Run a specific file interactively
include("notebooks/disaggregate_cube_minimal.jl")
```

There is no formal test suite. Development/exploration code lives in `src/function_dev.jl` and `notebooks/`.

## Architecture

The package is undergoing a major refactor on the `version_2.0_breaking` branch. The new architecture consolidates functionality into three utility files under `src/`, while the older modular structure in `src/datacube/` and `src/general/` is being phased out.

### Module Entry Point

`src/ItsLive.jl` exports:
- `disaggregate_cube` — **exported but not yet implemented**; cube-level disaggregation is currently done manually via `ItsLive.disaggregate` + `xmap`/`XOutput` (see `notebooks/disaggregate_cube_minimal.jl`)
- `dtbias_filter` — two-stage statistical filter (Spearman pre-screen + Mann-Whitney U test) for locking/skipping artifacts
- `sensor_bias_filter` — cross-sensor bias correction using grouped Mann-Whitney U tests
- `interval_average` — re-exported from TemporalDisaggregations; computes interval-averaged time series

### Internal API (accessed via `ItsLive.<name>`)

These are not exported but are the primary interface for data access and disaggregation:

| Function | Purpose |
|----------|---------|
| `datacube_catalog(catalog_geojson=datacube_catalog_path)` | Fetch GeoDataFrame of datacube polygons; URL or local path; defaults to S3 catalog |
| `datacube_path(lat, lon; catalog=datacube_catalog_path)` | Find zarr URL for a lat/lon point via spatial intersection; returns `missing` if none found |
| `datacube_load(path)` | Load zarr into RasterStack (lazy, via Zarr + YAXArrays) |
| `latlon_to_xy(lat, lon, rs)` | Convert lat/lon to projected x/y using Proj.jl (reads EPSG from `rs[:mapping]` metadata) |
| `disaggregate(vx, vy, t1, t2, sensor; ...)` | Per-pixel: apply filters + fit spline/GP/sinusoid; returns `DimStack` or `nothing` |
| `sensorgroup(sensor_labels)` | Map sensor strings → group IDs (Int16); always emits `@warn` about missing mission field |
| `plot_disaggregate(agg, t1, t2, disagg; ...)` | Dual-panel visualization; `disagg` can be a single `DimStack` or a vector of them |
| `to_datetime(ar; yeardecimal=false)` | Convert CFTime DimArray → Julia DateTime DimArray via `CFTime.timedecode` |
| `buffer_unidirectional(line, buffer; dims=1)` | Create offset polygon from LineString; `dims=2` buffers in the y-direction (used for velocity uncertainty corridor) |

Default catalog URL constant: `ItsLive.datacube_catalog_path = "https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json"`

### Source File Responsibilities

| File | Purpose |
|------|---------|
| `src/utilities_datacube.jl` | Core data access, `dtbias_filter`, `sensor_bias_filter`, `sensorgroup`, per-pixel `disaggregate` |
| `src/utilities_plot.jl` | `plot_disaggregate()` |
| `src/utilities.jl` | `to_datetime()`, `buffer_unidirectional()` |

### Legacy Source Files (src/datacube/, src/general/)

Not included in the module. Contain v1.x API: `getvar`, `lsqfit_annual`, `lsqfit_interp`, `sensorfilter`, `plotvar`, `plotbysensor`. Consult when reimplementing or porting.

### Key Data Flow

1. **Catalog lookup**: `datacube_catalog()` → GeoDataFrame → `datacube_path(lat, lon)` → zarr URL
2. **Data loading**: `datacube_load(path)` → `RasterStack`; or `YAXArrays.open_dataset(Zarr.zopen(path))` for `xmap`
3. **Coordinate transform**: `latlon_to_xy(lat, lon, rasterstack)` → projected x/y
4. **Per-pixel disaggregation**: `ItsLive.disaggregate(vx, vy, t1, t2, sensor; ...)` → `DimStack` with `:signal` and `:std`, or `nothing`
5. **Cube-level disaggregation**: `xmap(f, vx_cube ⊘ :mid_date, vy_cube ⊘ :mid_date, output=XOutput(...))` — apply `disaggregate` across all spatial pixels (see `notebooks/disaggregate_cube_minimal.jl`)

### Internal disaggregate Pipeline

`disaggregate` wraps two filters then calls `TemporalDisaggregations.disaggregate` twice:
0. Coarse validity screen (`abs(vx/vy) < max_abs_velocity=20000`, dates within output range); returns `nothing` if `< min_valid_observations=10` pass
1. `sensor_bias_filter` + `dtbias_filter` → combined keep mask; returns `nothing` if `< min_valid_observations` survive
2. First fit (with `loss_norm=:L1`) → builds a velocity-uncertainty corridor via `buffer_unidirectional(..., dims=2)`
3. Re-include filtered long-dt obs that intersect the fit corridor
4. Second fit on the expanded mask → final result (`DimStack` with `:signal` and `:std`)

### Cube-Level Disaggregation Pattern (xmap)

The canonical pattern (from `notebooks/disaggregate_cube_minimal.jl`) wraps `disaggregate` in a closure that captures shared metadata, and uses `Base.Fix` to partially apply the positional arguments `t1`, `t2`, `sensor`:

```julia
# Capture shared metadata (closed-over variables, computed once)
t1 = collect(ItsLive.to_datetime(ds["acquisition_date_img1"]))
t2 = collect(ItsLive.to_datetime(ds["acquisition_date_img2"]))
sensor = collect(ds["satellite_img1"])

output_start = Date(2014, 1, 1)
output_end   = Date(2025, 1, 1)
output_period = Week(1)

function f2(vx_col, vy_col)
    method = Spline(smoothness=1e-3, tension=10)
    # Base.Fix binds positional args 5, 4, 3 (sensor, t2, t1) onto disaggregate
    f = Base.Fix{3}(Base.Fix{4}(Base.Fix{5}(ItsLive.disaggregate, sensor), t2), t1)
    result = try
        f(vx_col, vy_col; method, output_start, output_end, output_period)
    catch e
        e isa LinearAlgebra.SingularException || rethrow()
        nothing
    end
    n = length(output_start:output_period:output_end)
    isnothing(result) ? fill(NaN, n) : result.signal.data
end

out0 = XOutput(Ti(output_start:output_period:output_end); outtype=Float64)
r = xmap(f2, vx0 ⊘ :mid_date, vy0 ⊘ :mid_date, output=out0, inplace=false)

# Write output to zarr
gen_cube = compute_to_zarr(Dataset(layer=r), "output.zarr", overwrite=true, max_cache=1e9)
```

Key points:
- `⊘ :mid_date` selects the `mid_date` dimension to map over (not the spatial dimensions)
- `Base.Fix{N}` binds the Nth positional argument; chained to bind args 3 (t1), 4 (t2), 5 (sensor)
- `compute_to_zarr` materializes the lazy `xmap` result to disk

### Sensor Group IDs

`sensorgroup()` maps sensor label strings to integer group IDs used by `sensor_bias_filter`. Sentinel-2 is always group 1 (the reference group):

| ID | Sensor |
|----|--------|
| 1 | Sentinel-2 (2A, 2B) — **reference group** |
| 2 | Landsat 8/9 |
| 3 | Sentinel-1 (1A, 1B) |
| 4 | Landsat 7 |
| 5 | Landsat 4/5 |

### External Dependencies of Note

- **Zarr + YAXArrays**: cloud zarr access (S3 without download); `xmap`/`XOutput`/`compute_to_zarr` for multi-pixel mapped computation
- **GeoDataFrames + GeoJSON + GeometryOps**: spatial catalog queries
- **Proj**: CRS transformation (lat/lon ↔ projected datacube coordinates)
- **TemporalDisaggregations**: core disaggregation algorithms (GP, Spline, Sinusoid); `interval_average` re-exported from here
- **DateFormats**: `yeardecimal()` for converting `DateTime` → decimal year (used throughout)
- **CairoMakie**: all visualization
- **CFTime**: handling non-standard calendar types in ITS_LIVE zarr time axes (DateTimeNoLeap, DateTime360Day, etc.)

### Known Issues / Active Development

- `src/junk.jl` documents a GeoDataFrames URL-reading issue; workaround in `datacube_catalog`: fetch via HTTP then pass body bytes to `GDF.read`
- `sensorgroup` emits a `@warn` about missing mission field — known TODO
- `disaggregate_cube` is exported but not yet defined; the intended interface is TBD
- `src/function_dev.jl` contains in-progress development code, not part of the module
- The catalog functions (`src/catalog.jl`, `src/datacube/catalog.jl`) were deleted; catalog logic is now in `utilities_datacube.jl`
