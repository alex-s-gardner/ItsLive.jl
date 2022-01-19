# the first time this is run you need to add packages
# using Pkg
# Pkg.add("Zarr")
# Pkg.add("ESDL")
# Pkg.add("ArchGDAL")
# Pkg.add("AWS")
# Pkg.add("YAXArrays")
# Pkg.add("Statistics")
# Pkg.add("Dates")
# Pkg.add("Plots")
# Pkg.add("IterTools")
# Pkg.add("Distributed")
# Pkg.add("JuliennedArrays")
# Pkg.add(PackageSpec(url="https://github.com/esa-esdl/ESDL.jl"))

# load in libraries
using AWS, Zarr, Statistics, Dates, Plots, IterTools, Distributed, JuliennedArrays, ESDL

# set up aws configuration
AWS.global_aws_config(AWSConfig(creds=nothing, region = "us-west-2"))

# path to datacube
path = Array{String,1}()
#push!(path,"http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")
push!(path,"s3://its-live-data/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")

function threaded_read(xin)
    xout = similar(xin)
    Threads.@threads for i in map(i->i.indices,Zarr.DiskArrays.eachchunk(xin))
        xout[i...] = xin[i...]
    end
    xout
end

i = 1;

# map into datacube
z = zopen(path[i], consolidated=true)

# map to specific variables
vx = z["vx"]
vy = z["vy"]
dt = z["date_dt"]

# Loop for each vertical column of chubked arrays
cs = Zarr.DiskArrays.eachchunk(vx).chunksize[1:2]
sz = size(vx,[1,2])
xout = zeros(sz)
stpi = cld(cs[1],1)
stpj = cld(cs[2],1)
nanmed(x) = median(skipmissing(x))

# read in image-pair time seperation 
dt = dt[:]
dt_max = [16 32 64 128 256 (2^15 - 1)];

function dtmaxfind(x::Base.ReshapedArray{Union{Missing, Int16}, 3, Array{Union{Missing, Int16}, 3}, Tuple{}}, dt::Vector{Float64}, dt_max::Matrix{Int64})

    # define a median absolut difference function
    function medmad(x::Array{Union{Missing, Int16}}) 
        if all(ismissing, x)
            medx = 0
            madx = 0
        else
            medx = round(Int,median(skipmissing(x)));
            madx = round(Int,median(skipmissing(abs.(x .- medx))));
        end
        return [medx, madx]
    end

    # find last valid bin or return zero
    function findlastorzero(x::BitVector)
        if any(x)
            x = findlast(x)
        else
            x = 0
        end
        return x
    end

    dt = round.(Int, dt)

    # sort index
    dt_sortperm = sortperm(dt);
    dt = dt[dt_sortperm]

    # find the last sorted value that is < 
    bin_ind = searchsortedlast.(Ref(dt), dt_max) 
    bin_ind = [1 bin_ind]

    x = x[:,:,dt_sortperm]

    binMad = zeros(Int16, size(x,1), size(x,2), length(dt_max))
    binMed = zeros(Int16, size(x,1), size(x,2), length(dt_max))

    size(binMad)
    for i in range(1, length(bin_ind)-1)
            foo = mapslices(medmad,x[:,:,bin_ind[i]:bin_ind[i+1]],dims=3);
            binMed[:,:,i] = foo[:,:,1]
            binMad[:,:,i] = foo[:,:,2]
    end

    # dtbin_mad_thresh: used to determine in dt means are significantly different
    dtbin_mad_thresh = 2; 

    # check if populations overlap (use first, smallest dt, bin as reference)
    minBound = binMed - binMad * dtbin_mad_thresh * 1.4826;
    maxBound = binMed + binMad * dtbin_mad_thresh * 1.4826;

    exclude = (minBound .> maxBound[:,:,1]) .| (maxBound .< minBound[1]);

    maxdt = mapslices(findlastorzero, exclude, dims = 3)
    noMax = maxdt .== 0
    maxdt[.~noMax] = dt_max[maxdt[.~noMax]]
    maxdt[noMax] .= (2^15 - 1)

    return maxdt
end

maxdt = zeros(Int16,sz)
Threads.@threads for i in range(1, step=stpi, length=cld(sz[1],stpj))
# for i in range(1, step=stpi, length=cld(sz[1],stpj))  
    j = for j in range(1, step=stpi, length=cld(sz[2],stpj))
        println([i,j])
 
        # Threaded read hangs - I think this has to do with the choice zarr chunk compression  
        # v2 = view(v, i:min(sz[1], i*stpi), j:min(sz[2], j*stpj), :)
       
        # threaded load into memory
        # @time v2 = threaded_read(v2)
        
        vxin = vx[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :]
    
        maxdt_vx = dtmaxfind(vxin, dt, dt_max)

        vyin = vy[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :]
        maxdt_vy = dtmaxfind(vyin, dt, dt_max)

        foo = cat(maxdt_vx, maxdt_vy,dims = 3);
        maxdt[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1)] = minimum(foo,dims=3)
    end
end

gr()
heatmap(1:size(maxdt,1),
    1:size(maxdt,2), maxdt,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="Median Velocity")