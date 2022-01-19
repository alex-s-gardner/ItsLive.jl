
# the first time this is run you need to add packages
#] add zarr
#] add AWS
#] add GLMakie

# load in libraries
using Zarr, AWS

# set up aws configuration
AWS.global_aws_config(AWSConfig(creds=nothing, region = "us-west-2"))

# path to datacube
path = Array{String,1}()
path = push!(path,"s3://its-live-data.jpl.nasa.gov/datacubes/v1/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")
path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_middate_11117.zarr")
path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_xy_10.zarr")
path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_xy_50_rechunker.zarr")
path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_xy100_with_rechunker.zarr")
path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_xy_10_create.zarr")
#path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_xy2_with_rechunker.zarr")
#path = push!(path,"s3://its-live-data.jpl.nasa.gov/test_datacube/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_xy1_with_rechunker.zarr")
#push!(path,"http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")

function threaded_read(xin)
    xout = similar(xin)
    Threads.@threads for i in map(i->i.indices,Zarr.DiskArrays.eachchunk(xin))
        xout[i...] = xin[i...]
    end
    xout
end

println("---- v[1, 1, :]-----")
for i = 1:size(path,1)
    t = time();
    # map into datacube
    z = zopen(path[i], consolidated=true)

    # map to specific variables
    v = z["v"]
    x = z["date_center"]

    # seperate file name and print
    foo = Base.Filesystem.splitdir(path[i])
    println(foo[2])

    # time read
    # read data using multiple threads
    v2 = view(v,1,1,:)

    a = threaded_read(v2)
    a = a == missing
    println(round(time()-t, digits = 1), " sec")
    println("--------------------------")
end

println("---- v[1:10, 1:10, :]-----")
for i = 1:size(path,1)
    t = time();
    # map into datacube
    z = zopen(path[i], consolidated=true)

    # map to specific variables
    v = z["v"]
    x = z["date_center"]



    # seperate file name and print
    foo = Base.Filesystem.splitdir(path[i])
    println(foo[2])

    # time read
    # read data using multiple threads
    v2 = view(v,1:10,1:10,:)
    a = threaded_read(v2)
    a = a == missing
    println(round(time()-t, digits = 1), " sec")
    println("--------------------------")
end


println("---- v[1:50, 1:50, :]-----")
for i = 1:size(path,1)
    t = time();
    # map into datacube
    z = zopen(path[i], consolidated=true)

    # map to specific variables
    v = z["v"]
    x = z["date_center"]



    # seperate file name and print
    foo = Base.Filesystem.splitdir(path[i])
    println(foo[2])

    # time read
    # read data using multiple threads
    v2 = view(v,1:50,1:50,:)
    a = threaded_read(v2)
    a = a == missing
    println(round(time()-t, digits = 1), " sec")
    println("--------------------------")
end


# threaded_read(v2)
#
# # plot data from
# @time foo = scatter(x,v[1,1,:]) #takes ~2 min
#
#
#
# @time b = v[1:100,1:100,:]; # takes 1.7 min
#
# image(a)
#
# idx = a .=== missing;
# x = x[:];
# a = Float32.(a[.!idx])
# x = Float32.(x[.!idx])
#
# spl = fit(SmoothingSpline, Float64.(x), Float64.(a), Float64.(250.0)) # λ=250.0
# Ypred = predict(spl) # fitted vector
# scene = scatter(x, a, color = :blue),scatter!(scene.axis,x,Ypred, color = :red)
#
#
# using SmoothingSplines
# using Gadfly
# using Zarr, AWS, GLMakie
# using Makie
# using SmoothingSplines
# using Gadfly
