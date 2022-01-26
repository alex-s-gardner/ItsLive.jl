"""
    datacube_map(catalogdf_row])

this function returns the zarr object and row/column for the zarr datacubes that intersects the provided latitude and longituded [decimal degrees].

use datacube_catalog.jl to generate the DataFrame catalog of the ITS_LIVE zarr datacubes and use datacube_intersect.jl to deteremine intersecting datacube

using Zarr

# Example
```julia
julia> datacube_map(69.1,-49.4, catalogdf)
```

# Arguments
   - `lat::Number`: latitude between -90 and 90 degrees
   - `lon::Number`: latitude between -180 and 180 degrees
   - `catalogdf::DataFrame`: DataFrame catalog of the ITS_LIVE zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasedena, California
January 25, 2022
"""