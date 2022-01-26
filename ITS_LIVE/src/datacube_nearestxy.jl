"""
datacube_nearestxy(lat,lon,dc)

this function returns the x/y indicies into a ZarrGroup for the points nearest the provided lat, lon locations

using Proj4

# Example no inputs
```julia
julia> datacube_catalog()
```

```julia
julia> datacube_catalog(catalog_geojson = "path/to/catalog.json")
```

# Arguments
   - `catalog_geojson::String`: path to geojson catalog of ITS_LIVE datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasedena, California
January 25, 2022
"""

function datacube_nearestxy(lat, lon, dc)
  # check that lat is within range
  if any((lat .< -90) .| (lat .> 90))
    error("lat = $lat, not in range [-90 to 90]")
  end

  # check that lon is within range
  if any((lon .<-180) .| (lon .> 180))
    error("lon = $lon, not in range [-180 to 180]")
  end

  function nearest_within(xPt, x)
    xWithin = abs(x[2] - x[1])
    xdist = abs.(x .- xPt)
    xdistmin,  xind0 = findmin(xdist)
  
    if (xdistmin > xWithin)
      xind0 = nothing
    end
  
    return xind0
  end
  
  # convert lat and lon into projected datacube coordinates
  trans = Proj4.Transformation("EPSG:4326", "EPSG:" * dc.attrs["projection"])
  x, y= trans.(eachrow(hcat(lat, lon)))

  # find the nearest x/y location with 1 full gridcell distance
  xind = nearest_within.(x, Ref(dc["x"][:]))
  yind = nearest_within.(y, Ref(dc["y"][:]))

  return xind, yind
end
