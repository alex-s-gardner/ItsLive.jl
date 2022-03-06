module ITS_LIVE

import DataFrames
import ArchGDAL
import Proj4
using AWS
using Zarr
using NamedArrays
using Dates
using Statistics
using DateFormats
using BSplineKit
using FastRunningMedian
using Polynomials
using Plots

include("datacube/catalog.jl")
include("datacube/intersect.jl")
include("datacube/nearestxy.jl")
include("datacube/getvar.jl")
include("datacube/dtfilter.jl")
include("datacube/vxvyfilter.jl")
include("datacube/decimalyear.jl")
include("datacube/lsqfit.jl")
include("datacube/lsqfit_interp.jl")
include("datacube/design_matrix.jl")
include("datacube/annual_matrix.jl")

end # module