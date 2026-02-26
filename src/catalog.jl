using GeoDataFrames, GeoJSON, HTTP

"""
    catalog([catalog_geojson::String])

Return a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes.
Loading is done via GeoDataFrames.jl; when `catalog_geojson` is a URL, it is
fetched with HTTP and then read from a temporary file.

# Examples

```julia
julia> catalog()
```

```julia
julia> catalog("https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json")
```

```julia
julia> catalog("path/to/catalog.json")
```

# Arguments

- `catalog_geojson::String`: URL or local path to the GeoJSON catalog of ITS_LIVE datacubes.

# Author

Alex S. Gardner, JPL, Caltech.
"""
function catalog(catalog_geojson::String = "https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json")
    if startswith(catalog_geojson, "http")
        response = HTTP.get(catalog_geojson)
        catalogdf = mktemp() do path, io
            write(io, response.body)
            close(io)
            GeoDataFrames.read(path)
        end
    else
        catalogdf = GeoDataFrames.read(catalog_geojson)
    end
    return catalogdf
end
