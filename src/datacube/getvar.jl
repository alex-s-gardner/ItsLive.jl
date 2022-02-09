"""
    getvar(lat,lon, varnames, catalogdf])

this function returns a named m x n matrix of vectors with m = length(lat) rows and n = length(varnames)+2(for lat and lon) columns for the points nearest the lat/lon location from ITS_LIVE Zarr datacubes

use catalog.jl to generate the DataFrame catalog of the ITS_LIVE zarr datacubes

using DataFrames

# Example
```julia
julia> getvar(69.1,-49.4, ["mid_date", "v"], catalogdf)
```

# Arguments
   - `lat::Union{Vector,Number}`: latitude between -90 and 90 degrees
   - `lon::Union{Vector,Number}`: latitude between -180 and 180 degrees
   - `varnames::Unions{String, Vector{String}}`: name of variables to extract from Zarr DataFrame
   - `catalogdf::DataFrame`: DataFrame catalog of the ITS_LIVE zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 4, 2022
"""
function getvar(lat::Union{Vector,Number},lon::Union{Vector,Number}, varnames::Union{String, Vector{String}}, catalogdf)

# check that lat is within range
if any(lat .<-90 .|| lat .> 90)
    error("lat = $lat, not in range [-90 to 90]")
end

# check that lon is within range
if any(lon .<-180 .|| lon .> 180)
    error("lon = $lon, not in range [-180 to 180]")
end

# check that catalog is a dataframe
if ~(catalogdf isa DataFrames.DataFrame)
    error("provided catalog is not a DataFrame, use catalog.jl to generate a DataFrame")
end

if(typeof(varnames) == String)
    varnames = [varnames]
end


# find the DataFrame rows of the datacube that intersect a series of lat/lon points
rows = ITS_LIVE.intersect.(lat,lon, Ref(catalogdf))

# find all unique datacubes (mutiple points can intersect the same cube)
urows = unique(rows)

rind = zeros(size(lat))
cind = zeros(size(lat))
vout = Vector{Vector{Union{Missing, Any}}}()
ind = Vector{Int64}()
vind = Vector{Int64}()

# select variable to extract
# varnames = ["vx"]
for row in urows
    # check if row is "missing"
    if ismissing(row)
        ind0 = findall(ismissing.(rows))

        for j = 1:lastindex(varnames)
            for i = 1:lastindex(ind0)
                push!(vout, [missing])
                push!(ind, ind0[i])
                push!(vind, j)
            end
        end
        continue
    end

    # extract path to Zarr datacube file
    path2cube = catalogdf[row,"zarr_url"]

    # load Zarr group
    dc = Zarr.zopen(path2cube)

    # find lat/lon values that intersect datacube
    ind0 = findall(x -> x==row, skipmissing(rows))

    # find closest point
    c,r = ITS_LIVE.nearestxy(lat[ind0], lon[ind0], dc)
    rind[ind0] .= r
    cind[ind0] .= c
    
    # extract timeseries from datacube

    # loop for each r and c
    for j = 1:lastindex(varnames)

        if ndims(dc[varnames[j]]) == 1
            foo = dc[varnames[j]][:]
            for i = 1:lastindex(r)
                push!(vout, foo)
                push!(ind, ind0[i])
                push!(vind, j)
            end
        elseif ndims(dc[varnames[j]]) == 3
            for i = 1:lastindex(r)
                push!(vout, dc[varnames[j]][r[i], c[i], :])
                push!(ind, ind0[i])
                push!(vind, j)
            end
        else
            error([varnames[j] " is not a 1D or 3D variable"])
        end
    end
end 

# arrange for rows of location (lat/lon) and columns of variables

# sort vector by variables
i = sortperm(vind)
ind = ind[i]
vout = vout[i]

# reshape putting unique variables into columns
vout = reshape(vout, (div(length(vout),length(varnames)), length(varnames)))

# sort rows to match input lat and lon
i = sortperm(ind[1:div(length(vout),length(varnames))])
vout = vout[i,:]

# add lat/lon to matrix
vout = hcat(lat,lon,vout)

# add naming to matrix
vout = NamedArrays.NamedArray(vout)
NamedArrays.setnames!(vout, vcat("lat","lon",varnames), 2)

return vout 
end