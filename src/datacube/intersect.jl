"""
    rownumber = intersect(lat,lon, catalogdf])

return the `rownumber` of the the DataFrame catalog (`catalogdf`) of the ITS_LIVE zarr datacubes 
that intersects the provided `lat`itude and `lon`gituded [decimal degrees].

use catalog.jl to generate the DataFrame catalog of the ITS_LIVE zarr datacubes

using ArchGDAL, DataFrames

# Example
```julia
julia> intersect(69.1,-49.4, catalogdf)
```

# Arguments
   - `lat::Number`: latitude between -90 and 90 degrees
   - `lon::Number`: latitude between -180 and 180 degrees
   - `catalogdf::DataFrame`: DataFrame catalog of the ITS_LIVE zarr datacubes
"""
function intersect(lat::Number,lon::Number, catalogdf::DataFrame)
# set up aws configuration

    # check that lat is within range
    if lat <-90 || lat > 90
        error("lat = $lat, not in range [-90 to 90]")
    end

    # check that lon is within range
    if lon <-180 || lon > 180
        error("lon = $lon, not in range [-180 to 180]")
    end

    # define a point
    point = ArchGDAL.createpoint(lon, lat)

    # find intersecting datacube polygon
    found = false
    for row in eachrow(catalogdf)
        if ArchGDAL.contains(row[1], point)
            found = true
            return DataFrames.rownumber(row)
            break
        end
    end

    # check if a datacube polygon was found
    if ~found
        @warn "could not find intersecting datacube for point[lat = $lat, lon = $lon]"
        return  missing
    end    
end