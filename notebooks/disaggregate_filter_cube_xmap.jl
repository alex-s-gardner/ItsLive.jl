# # disaggregate_filter_cube_xmap.jl
#
# Cube-level disaggregation using YAXArrays xmap pattern.
# Demonstrates:
#   • Cloud-native zarr access with lazy evaluation
#   • Chunk-aware processing via xmap
#   • Direct-to-disk output via savecube
#   • Sensor-bias and dt-bias filtering per pixel
#   • Many-to-many pattern: single xmap call outputs both vx and vy separately (no redundant computation!)
#
# Contrast with disaggregate_filter_cube_example.jl which uses manual iteration.
#
# USAGE:
#   # Run with default worker count (7):
#   julia --project=. notebooks/disaggregate_filter_cube_xmap.jl
#
#   # To adjust worker count, edit n_workers in the script
#   # Note: Uses Distributed.jl (multi-process) for parallelism, not threading
#
# Author: Alex S. Gardner, JPL / Caltech.

# ── 1. Dependencies ─────────────────────────────────────────────────────────

# Set up distributed computing FIRST (before loading heavy packages)
using Distributed
if nworkers() == 1
    n_workers = 7  # Adjust based on your CPU cores
    println("Adding $n_workers worker processes...")
    addprocs(n_workers)
    println("✓ Total processes: ", nprocs(), " (1 main + ", nworkers(), " workers)")
end

# Load packages on ALL processes (main + workers)
@everywhere begin
    using ItsLive
    using TemporalDisaggregations
    using YAXArrays
    using Zarr
    using DimensionalData
    using Dates
    using Statistics
    using DateFormats: yeardecimal
    using KernelFunctions
    using LinearAlgebra

    # For faster BLAS on macOS:
    using AppleAccelerate

    # For faster BLAS on Intel Linux:
    # using MKL

    # Zarr.MaxLengthStrings.length() returns `nothing` for all-zero (empty) fill-value strings,
    # which propagates into ncodeunits and breaks JSON serialization via `1:nothing`.
    # Fix: return 0 instead of nothing for empty MaxLengthStrings.
    Base.ncodeunits(s::Zarr.MaxLengthStrings.MaxLengthString) = something(findlast(!iszero, s.data), 0)
end

# CairoMakie only needed on main process for plotting
using CairoMakie

# ── 2. Disaggregation configuration (on all workers) ────────────────────────
@everywhere begin
    output_period = Month(1)
    loss_norm = :L1
    sigma_buffer = 2
    time_buffer = Month(1)
    observation_error_minimum = 20 # m/yr
    verbose = false

    method = Spline(smoothness=1e-3, tension=1)
end

# ── 3. Locate and load datacube ─────────────────────────────────────────────
begin
    lat, lon = 60.0626, -139.3193   # Hubbard Glacier, Alaska
    zarr_path = ItsLive.datacube_path(lat, lon)
    local_dir = "/Users/gardnera/data/its-live-data/datacubes/v02"
    zarr_path = joinpath(local_dir, last(splitpath(zarr_path)))
    println("Zarr path: ", zarr_path)
end

# Load as YAXArrays Dataset (required for xmap)
begin
    ds = YAXArrays.open_dataset(Zarr.zopen(zarr_path))
    println("Dataset loaded with dimensions: ", keys(ds.axes))
    println("Available variables: ", keys(ds.cubes))

    # Report parallel configuration
    println("\n⚡ Parallel configuration:")
    println("   Workers: ", nworkers())
    println("   Total processes: ", nprocs(), " (1 main + ", nworkers(), " workers)")
    println("   Expected speedup: ~", min(nworkers(), 4), "× for large datasets")
end

# ── 4. Extract and prepare metadata (distribute to all workers) ────────────
begin
    # Extract time series metadata on main process
    t1_main = collect(ItsLive.to_datetime(ds["acquisition_date_img1"]))
    t2_main = collect(ItsLive.to_datetime(ds["acquisition_date_img2"]))
    sensor_main = String.(collect(ds["satellite_img1"]))

    # Create error arrays (placeholder - using uniform errors as in original)
    vx_err_main = ones(size(t1_main))
    vy_err_main = ones(size(t1_main))

    # Set minimum error
    vx_err_main[vx_err_main .< observation_error_minimum] .= observation_error_minimum
    vy_err_main[vy_err_main .< observation_error_minimum] .= observation_error_minimum

    # Sort data by t1 if needed
    if !issorted(t1_main)
        order = sortperm(t1_main)
        t1_main = t1_main[order]
        t2_main = t2_main[order]
        vx_err_main = vx_err_main[order]
        vy_err_main = vy_err_main[order]
        sensor_main = sensor_main[order]
        println("Data sorted by acquisition date")
    end

    # Define output time range
    output_start = Date(2014, 1, 1) + time_buffer
    output_end = Date(maximum(skipmissing(t2_main))) - time_buffer
    time_output = output_start:output_period:output_end
    n_time = length(time_output)

    println("Output time range: $(output_start) to $(output_end)")
    println("Number of output time steps: $(n_time)")
    println("Number of input observations: $(length(t1_main))")

    # Distribute metadata to ALL workers
    println("Distributing metadata to all workers...")
    @everywhere t1 = $t1_main
    @everywhere t2 = $t2_main
    @everywhere sensor = $sensor_main
    @everywhere vx_err = $vx_err_main
    @everywhere vy_err = $vy_err_main
    @everywhere output_start = $output_start
    @everywhere output_end = $output_end
    @everywhere output_period = $output_period
    println("✓ Metadata distributed to ", nworkers(), " workers")
end

# ── 5. Define xmap wrapper function (on all workers) ───────────────────────
@everywhere begin
    # Wrapper function using Many InDims to many OutDims pattern
    # Outputs FIRST (vx_out, vy_out), then inputs (vx_col, vy_col)
    # Mutates outputs in-place and returns nothing
    function disagg_wrapper_separate!(vx_out, vy_out, vx_col, vy_col)
        # Handle missing data
        if all(ismissing.(vx_col)) || all(ismissing.(vy_col))
            vx_out .= NaN
            vy_out .= NaN
            return nothing
        end

        # Call disaggregate with all closed-over variables
        # These variables (method, vx_err, vy_err, t1, t2, sensor, etc.) are now
        # available on all workers thanks to @everywhere blocks above
        result = try
            ItsLive.disaggregate(
                method, vx_col, vy_col, vx_err, vy_err, t1, t2, sensor;
                output_start, output_end, output_period,
                loss_norm, sigma_buffer, time_buffer, verbose)
        catch e
            # Handle singular matrix exceptions (insufficient data)
            if e isa LinearAlgebra.SingularException
                vx_out .= NaN
                vy_out .= NaN
                return nothing
            else
                rethrow()
            end
        end

        vx_fit, vy_fit, valid_obs = result

        # Mutate outputs in-place
        if isnothing(vx_fit)
            vx_out .= NaN
            vy_out .= NaN
        else
            vx_out .= vx_fit.signal.data
            vy_out .= vy_fit.signal.data
        end

        return nothing  # Required for many-to-many pattern
    end
end

println("✓ Wrapper function defined on all workers ")

# ── 6. Apply xmap (single call for both vx and vy) ─────────────────────────
begin
    # Subset for faster testing (comment out to process full cube)
    TEST_SUBSET = true

    # Get velocity cubes and optionally subset
    if TEST_SUBSET
        println("\n⚠ TEST MODE: Processing subset (10x10 = 100 pixels)")
        # Use trailing subset to avoid missing chunks and faster testing
        vx_cube = ds.vx[end-99:end, end-99:end, :] ⊘ :mid_date
        vy_cube = ds.vy[end-99:end, end-99:end, :] ⊘ :mid_date
    else
        vx_cube = ds.vx ⊘ :mid_date
        vy_cube = ds.vy ⊘ :mid_date
    end

    println("Input cube dimensions: ", dims(ds["vx"]))
    println("Mapping over dimension: mid_date")

    # Define separate output specifications for vx and vy
    # xmap will automatically add spatial dimensions from input cubes
    out_vx = YAXArrays.XOutput(Ti(time_output); outtype=Float64)
    out_vy = YAXArrays.XOutput(Ti(time_output); outtype=Float64)

    # Apply xmap ONCE for both components using many-to-many pattern
    println("\nApplying xmap with separate vx/vy outputs...")
    println("  Input vx shape: ", size(vx_cube))
    println("  Input vy shape: ", size(vy_cube))

    (result_vx, result_vy) = xmap(
        disagg_wrapper_separate!,
        vx_cube,
        vy_cube,
        output=(out_vx, out_vy),
        inplace=true  # Required for many-to-many pattern
    )

    println("✓ xmap lazy computation graph created")
    println("  Result vx dimensions: ", dims(result_vx))
    println("  Result vy dimensions: ", dims(result_vy))
    println("  Result vx shape: ", size(result_vx))
    println("  Result vy shape: ", size(result_vy))
    println("✓ Each pixel processed ONCE (no redundant computation!)")
end

# ── 7. Materialize to zarr ──────────────────────────────────────────────────
begin
    output_path = joinpath(@__DIR__, "output_disagg_xmap.zarr")

    println("\nWriting to: $(output_path)")
    println("Computing output...")
    dsout = Dataset(vx=result_vx, vy=result_vy)

    @time compute_to_zarr(dsout, output_path, overwrite=true, max_cache=5000)

    println("✓ computation complete!")
    println("Outputs written to: $(output_path)")

end

# ── 8. Verify output ───────────────────────────────────────────────────────
begin
    # Load output dataset
    ds_out = YAXArrays.open_dataset(Zarr.zopen(output_path))

    println("\nOutput verification:")
    println("  Variables: ", keys(ds_out.cubes))
    println("  vx shape: ", size(ds_out["vx"]))
    println("  vy shape: ", size(ds_out["vy"]))
    println("  Dimensions: ", keys(ds_out.axes))

    # Check for non-NaN values (sample first 10 time steps)
    vx_sample = collect(ds_out["vx"][Ti=1:min(10, n_time)])
    vy_sample = collect(ds_out["vy"][Ti=1:min(10, n_time)])

    n_valid_vx = count(!isnan, vx_sample)
    n_valid_vy = count(!isnan, vy_sample)
    n_total = length(vx_sample)

    println("\n  Valid vx pixels (first 10 times): $(n_valid_vx) / $(n_total) ($(round(100*n_valid_vx/n_total, digits=1))%)")
    println("  Valid vy pixels (first 10 times): $(n_valid_vy) / $(n_total) ($(round(100*n_valid_vy/n_total, digits=1))%)")

    if n_valid_vx > 0
        println("\n  vx stats: mean = $(round(mean(skipmissing(vx_sample)), digits=2)) m/yr")
        println("  vy stats: mean = $(round(mean(skipmissing(vy_sample)), digits=2)) m/yr")
    end
end

println("\n✓✓✓ Script complete! ✓✓✓")
println("\nAdvantages over manual iteration:")
println("  ✓ No pre-allocation in memory (streams to disk)")
println("  ✓ Chunk-aware processing (respects zarr chunk boundaries)")
println("  ✓ Scalable to arbitrarily large cubes")
println("  ✓ Distributed parallel processing across ", nworkers(), " workers")
println("  ✓ Single xmap call (no redundant computation!)")
println("  ✓ Separate vx/vy outputs (no concatenation/splitting needed!)")
println("\nTo load output:")
println("  ds = YAXArrays.open_dataset(Zarr.zopen(\"$(output_path)\"))")
println("  vx = ds[\"vx\"]")
println("  vy = ds[\"vy\"]")
