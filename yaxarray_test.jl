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
dt = z["dt"]

# Loop for each vertical column of chubked arrays
cs = Zarr.DiskArrays.eachchunk(v).chunksize[1:2]
sz = size(v,[1,2])
xout = zeros(sz)
stpi = cld(cs[1],1)
stpj = cld(cs[2],1)
nanmed(x) = median(skipmissing(x))

#Threads.@threads for i in range(1, step=cs[1], length=cld(sz[1],cs[1])-1)
for i in range(1, step=stpi, length=cld(sz[1],stpj))
    
    j = for j in range(1, step=stpi, length=cld(sz[2],stpj))
        println([i,j])
 
        #return
        #v2 = view(v, i:min(sz[1], i*stpi), j:min(sz[2], j*stpj), :)
       
        # load data into memory
        # @time v2 = threaded_read(v2)
        @time xin = v[i:min(sz[1], i*stpi), j:min(sz[2], j*stpj), :]

        #@time xin = threaded_read(view(v, i:min(sz[1], i*stpi), j:min(sz[2], j*stpj), :))

        @time xout[i:min(sz[1], i*stpi), j:min(sz[2], j*stpj)] = map(nanmed, Slices(xin, 3))
    end
end

return

gr()
heatmap(1:size(xout,1),
    1:size(xout,2), xout,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="Median Velocity")
#
path = ["s3://its-live-data/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr"];
ds = open_dataset(path[1])
#using Distributed
#addprocs(4)
#@everywhere using Statistics
#@everywhere using ESDL

# a = ds.v
indims = InDims(getAxis("mid_date",ds.v))
outdims = OutDims(getAxis("mid_date",ds.v))
mapCube(mean, ds.v, indims=indims, outdims=outdims)
# @time mapslices(mean, ds.v.data, dims = "mid_date")

@time f = a[:,:,range(1,250)]

# function to conver between datetime64 and julia datetime
const days2ns = 86400000000000.56
npdt64_to_dt(t) = Dates.Nanosecond(round(t*days2ns)) + DateTime(1970)

# convert from decimal day to julia datetime
t = npdt64_to_dt.(ds.mid_date.values[:])

# take median of all data in a column in parallel
idx = t .> DateTime(2018)
ds.mid_date.values[idx]

@time  f = a[range(1,100),range(1,100),range(1,250)]
@time s.mean(f[:,:,5])

for i in eachindex(view(f, 1:100, 1:100,1))
    println(i)
end



function cubefilt(x, dt, dt_edge, dtbin_mad_thresh)
    # filter data cube

    # dtbin_mad_thresh: used to determine in dt means are significantly different
    # dtbin_mad_thresh = 2*1.4826; 

    # dt_edge: edges of dt bins
    # dt_edge = [0 16 32 64 128 256 inf];

    # check if populations overlab (use first, smallest dt, bin as reference)
    [m, n, ~] = size(x);

% initialize output and functions
invalid = false(size(x));
maxdt = nan([m,n]);
madFun = @(x) median(abs(x - median(x))); % replaced on chad's suggestion

% parfor version of loop
parfor i = 1:m
    % slice for parfor
    x1 = x(i,:,:);
    maxdt1 = maxdt(i,:);
    invalid1 = invalid(i,:,:);
    for j = 1:n
        x0 = squeeze(x1(1,j,:));
        if all(isnan(x0))
            continue
        end

        %% are means significantly different for various dt groupings?
        ind = discretize(dt(~isnan(x0)),dt_edge);
        [xm, xedg] = groupsummary(x0(~isnan(x0)),ind,'median');
        xmad = groupsummary(x0(~isnan(x0)),ind,madFun);
        
        % check if populations overlap (use first, smallest dt, bin as reference)
        minBound = xm - xmad*dtbin_mad_thresh* 1.4826;
        maxBound = xm + xmad*dtbin_mad_thresh* 1.4826;
        
        exclude = (minBound > maxBound(1)) | (maxBound < minBound(1));
        
        if any(exclude)
            maxdt1(1,j) = min(dt_edge(xedg(exclude)));
            invalid1(1,j,:) = dt > maxdt1(1,j);
        end
       
    end
    % place into larger array for parfor
    maxdt(i,:) = maxdt1;
    invalid(i,:,:) = invalid1;
end



