module ITS_LIVE

import DataFrames
import ArchGDAL
import Proj4
using AWS
using Zarr
using NamedArrays
using Statistics
using DateFormats
using BSplineKit
using FastRunningMedian
using Polynomials
using Plots
using NearestNeighbors

include("datacube/catalog.jl")
include("datacube/intersect.jl")
include("datacube/nearestxy.jl")
include("datacube/getvar.jl")
include("datacube/dtfilter.jl")
include("datacube/vxvyfilter.jl")
include("datacube/decimalyear.jl")
include("datacube/lsqfit.jl")
include("datacube/lsqfit_itslive.jl")
include("datacube/lsqfit_interp.jl")
include("datacube/design_matrix.jl")
include("datacube/annual_matrix.jl")
include("datacube/wlinearfit.jl")
include("datacube/running_mean.jl")
include("datacube/annual_magnitude.jl")
include("datacube/climatology_magnitude.jl")

end # module