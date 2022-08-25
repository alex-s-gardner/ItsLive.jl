module ItsLive

using DataFrames
using ArchGDAL
using Proj4
using AWS
using Zarr
using NamedArrays
using OrderedCollections
using Statistics
using DateFormats
using BSplineKit
using FastRunningMedian
using Polynomials
using Plots
using NearestNeighbors

include("general/binstats.jl")
include("general/running_mean.jl")
include("general/decimalyear.jl")

include("datacube/catalog.jl")
include("datacube/intersect.jl")
include("datacube/nearestxy.jl")
include("datacube/getvar.jl")
include("datacube/dtfilter.jl")
include("datacube/vxvyfilter.jl")
include("datacube/lsqfit_annual.jl")
include("datacube/lsqfit_interp.jl")
include("datacube/design_matrix.jl")
include("datacube/annual_matrix.jl")
include("datacube/wlinearfit.jl")
include("datacube/annual_magnitude.jl")
include("datacube/climatology_magnitude.jl")
include("datacube/sensorgroup.jl")
include("datacube/plotbysensor.jl")
include("datacube/sensorfilter.jl")
include("datacube/plotvar.jl")
include("datacube/save2h5.jl")

# precompile functions
precompile(plotvar,(NamedMatrix{Any, Matrix{Any}, Tuple{OrderedDict{String, Int64}, OrderedDict{String, Int64}}}, String, Number))
precompile(getvar, (Union{Vector,Number},Union{Vector,Number}, Union{String, Vector{String}}, DataFrame))
precompile(catalog, ())
precompile(catalog, (String))

end # module