"""
    datacube_intersect(lat,lon, catalogdf])

this function returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longituded [decimal degrees].

use datacube_catalog.jl to generate the DataFrame catalog of the ITS_LIVE zarr datacubes

using ArchGDAL, DataFrames

# Example
```julia
julia> datacube_intersect(69.1,-49.4, catalogdf)
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

function datacube_intersect(lat::Number,lon::Number, catalogdf)
# set up aws configuration

    # check that lat is within range
    if lat <-90 || lat > 90
        error("lat = $lat, not in range [-90 to 90]")
    end

    # check that lon is within range
    if lon <-180 || lon > 180
        error("lon = $lon, not in range [-180 to 180]")
    end

    # check that catalog is a dataframe
    if ~(catalogdf isa DataFrame)
        error("provided catalog is not a DataFrame, use datacube_catalog.jl to generate a DataFrame")
    end

    # define a point
    point = ArchGDAL.createpoint(lon, lat)

    # find intersecting datacube polygon
    found = false
    for row in eachrow(catalogdf)
        if ArchGDAL.contains(row[1], point)
            found = true
            return rownumber(row)
            break
        end
    end

    # check if a datacube polygon was found
    if ~found
        @warn "could not find intersecting datacube for point[lat = $lat, lon = $lon]"
        return  missing
    end    
end