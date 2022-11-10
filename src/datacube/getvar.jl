"""
    df = getvar(lat, lon, varnames, catalogdf; s3 = false)

return a DataFrame (`df`) with m = length(`lat`) rows 
and n = length(`varnames`)+2(for `lat` and `lon`) columns for the points nearest 
the `lat`/`lon` location from ITS_LIVE Zarr datacubes

use `catalog.jl` to generate the DataFrame catalog (`catalogdf`) of the ITS_LIVE zarr datacubes

using DataFrames Dates

# Example
```julia
julia> getvar(69.1,-49.4, ["mid_date", "v"], catalogdf)
```

# Arguments
   - `lat::Union{Vector,Number}`: latitude between -90 and 90 degrees
   - `lon::Union{Vector,Number}`: latitude between -180 and 180 degrees
   - `varnames::Unions{String, Vector{String}}`: name of variables to extract from Zarr DataFrame
   - `catalogdf::DataFrame`: DataFrame catalog of the ITS_LIVE zarr datacubes
   - `s3 =false`: keyword argument for using s3 vs https paths [optional]

"""
function getvar(lat::Union{Vector,Number},lon::Union{Vector,Number}, varnames::Union{String, Vector{String}}, catalogdf::DataFrame; s3 = false)
    
# check that lat is within range
if any(lat .<-90) || any(lat .> 90)
    error("lat = $lat, not in range [-90 to 90]")
end

# check that lon is within range
if any(lon .<-180) || any(lon .> 180)
    error("lon = $lon, not in range [-180 to 180]")
end

# check that lon are the same length
if length(lon) != length(lat)
    error("lon and lat are not the same length")
end

# check that catalog is a dataframe
if ~(catalogdf isa DataFrames.DataFrame)
    error("provided catalog is not a DataFrame, use catalog.jl to generate a DataFrame")
end

if(typeof(varnames) == String)
    varnames = [varnames]
end

if ~isa(lat, Array)
    lat = [lat]
    lon = [lon]
end

# find the DataFrame rows of the datacube that intersect a series of lat/lon points
rows = ItsLive.intersect.(lat,lon, Ref(catalogdf))

# extract urls
zarr_url = catalogdf[!,"zarr_url"]
if s3
    zarr_url = replace.(zarr_url, "http://its-live-data.s3.amazonaws.com" => "s3://its-live-data")
end

# find all unique datacubes (mutiple points can intersect the same datacube)
urows = unique(rows)

C = Vector{Vector{Vector{Union{Missing, Any}}}}()
ind = Vector{Vector{Int64}}()
vind = Vector{Vector{Int64}}()

for i in 1:Threads.nthreads()
    push!(C,Vector{Vector{Union{Missing, Any}}}())
    push!(ind,Vector{Int64}())
    push!(vind,Vector{Int64}())
end

# parallel loop over variables (not rows) is most efficient in nearly all cases
for row in urows
    # check if row is "missing"
    if ismissing(row)
        ind0 = findall(ismissing.(rows))

        for j in eachindex(varnames)
            for i in eachindex(ind0)
                push!(C[1], [missing])
                push!(ind[1], ind0[i])
                push!(vind[1], j)
            end
        end
        continue
    end

    # extract path to Zarr datacube file
    path2cube = zarr_url[row]

    # load Zarr group
    dc = Zarr.zopen(path2cube; fill_as_missing = true)

    # find lat/lon values that intersect datacube
    ind0 = findall(x -> x==row, skipmissing(rows))

    # find closest point 
    # NOTE: in Zarr cube "x" changes along rows (r) and "y" changes along columns (c)
    cartind = ItsLive.nearestxy(lat[ind0], lon[ind0], dc)
    
    ## print row and column
    # println(path2cube)
    # println("row = ", r, ", col = ", c)

    # extract timeseries from datacube
    
    # loop for each variable
    Threads.@threads for j in eachindex(varnames)
    #for j in eachindex(varnames)
        if ndims(dc[varnames[j]]) == 1
            foo = dc[varnames[j]][:]
            for i in eachindex(cartind)
                push!(C[Threads.threadid()], foo)
                push!(ind[Threads.threadid()], ind0[i])
                push!(vind[Threads.threadid()], j)
            end
        elseif ndims(dc[varnames[j]]) == 3
            push!(C[Threads.threadid()], dropdims(dc[varnames[j]][cartind,:]; dims=1))

            for i in eachindex(cartesianind)
                push!(ind[Threads.threadid()], ind0[i])
                push!(vind[Threads.threadid()], j)
            end
        else
            error([varnames[j] " is not a 1D or 3D variable"])
        end
    end
end

# concat idividual thread arrays
C = vcat(C...)
ind = vcat(ind...)
vind = vcat(vind...)

# arrange for rows of location (lat/lon) and columns of variables

# sort vector by variables
i = sortperm(vind)
ind = ind[i]
C = C[i]

# reshape putting unique variables into columns
C = reshape(C, (div(length(C),length(varnames)), length(varnames)))

# sort rows to match input lat and lon
i = sortperm(ind[1:div(length(C),length(varnames))])
C = C[i,:]

# add lat/lon to matrix
C = hcat(lat,lon,C)

# create a DataFrame
df = DataFrame()
vars = vcat("lat", "lon",varnames);

# find datetime variables and convert 
datevarnames = ["acquisition_date_img2", "acquisition_date_img1", "date_center", "mid_date"]

# define a function to convert from python datetime to Julia datetime

# create SecondMissing function to account for missing dates
SecondMissing(x::Missing) = missing; 
SecondMissing(x::Number) =  Dates.Second(x); 
npdt64_to_dt(t) =  SecondMissing.((round.(t.*86400))) .+ Dates.DateTime(1970) 

for j in findall(in(datevarnames),vars)
    C[:,j] = [npdt64_to_dt.(row[1]) for row = eachrow(C[:,j])]
end


# define a function to convert fixed type.. also convert Zarr.MaxLengthStrings.MaxLengthString
# to Stings as Zarr.MaxLengthStrings to play nice with other packages like Arrow.jl
function convertvec(x)
    if (x[1] isa Zarr.MaxLengthStrings.MaxLengthString{1, UInt32}) || 
        (x[1] isa Zarr.MaxLengthStrings.MaxLengthString{2, UInt32}) ||
        (x[1] isa Zarr.MaxLengthStrings.MaxLengthString{3, UInt32}) ||
        (x[1] isa Zarr.MaxLengthStrings.MaxLengthString{32, UInt32})

        if any(ismissing.(x))
            x = convert(Union{Missing,Vector{String}}, x) 
        end
            x = convert(Vector{String}, x)
    end
    else
        if any(ismissing.(x))
            x = convert(Union{Missing, Vector{typeof(x[1])}}, x)
        else
            x = convert(Vector{typeof(x[1])}, x)
        end
    end
end

# convert data types
for i in eachindex(vars)
    if any(vars[i] .== ["lat", "lon"])
        df[!,vars[i]] = convert(Vector{Float64},C[:,i])
    else
        df[!,vars[i]] =  [convertvec(row[1]) for row in eachrow(C[:,i])]; #convert.(Float64,C[:,i])
    end
end

return df
end