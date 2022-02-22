# test read speads and memory overhead for various datacube structures
using Zarr

n = 100;
r = rand(1:833, n);
c = rand(1:833, n);

print("row = ")
println(r)

print("col = ")
println(c)

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_t_25730_xy_10_test.zarr" 
#Cube sizes: (25730, 833, 833)
#Chunk sizes: (25730, 10, 10)
#Access v[1, 1]: 0.3554554230000804 seconds
#Access v[5, 5]: 0.32438522399991143 seconds
#Access v[400, 400]: 0.32991992699999173 seconds
#Total Objects: 68264
#Total Size: 12.9 GiB
#S3 PUT cost: $0.34
#Runtime: estimated at 13-18 hours as original file is rechunked from (250,100,100)
println(path2cube);
dc = Zarr.zopen(path2cube)
@time for i = 1:n; dc["v"][r[i],c[i],:]; end
# 15.436011 seconds (749.70 k allocations: 571.617 MiB, 0.46% gc time, 2.00% compilation time)

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_t_25730_xy_5_test.zarr" 
#Cube sizes: (25730, 833, 833)
#Chunk sizes: (25730, 5, 5)
#Access v[1, 1]: 0.3262468899999931 seconds
#Access v[5, 5]: 0.3098814530000027 seconds
#Access v[400, 400]: 0.3079671650000364 seconds
#Total Objects: 255738
#Total Size: 15.4 GiB
#S3 PUT price: $1.27
#Runtime: estimated at 30 hours as original file is rechunked from (250,100,100)
println(path2cube);
dc = Zarr.zopen(path2cube)
@time for i = 1:n; dc["v"][r[i],c[i],:]; end
# 6.979865 seconds (11.25 k allocations: 136.609 MiB, 0.34% gc time)

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_t_10K_xy_5_complete.zarr"
#Cube sizes: (25730, 833, 833)
#Chunk sizes: (10000, 5, 5)
#Access v[1, 1]: 0.3341476059999877 seconds
#Access v[5, 5]: 0.277660626999932 seconds
#Access v[400, 400]: 0.34147050599995055 seconds
#Total Objects: 757803
#Total Size: 15.4 GiB
#S3 PUT price: $3.78
#Runtime: 14 hours
println(path2cube);
dc = Zarr.zopen(path2cube)
@time for i = 1:n; dc["v"][r[i],c[i],:]; end

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_t_10K_xy_10_complete.zarr"
#Cube sizes: (25730, 833, 833)
#Chunk sizes: (10000, 10, 10)
#Access v[1, 1]: 0.5651649210012692 seconds
#Access v[5, 5]: 0.4418043139994552 seconds
#Access v[400, 400]: 0.4450376610002422 seconds
#Total Objects: 195402
#Total Size: 12.9 GiB
#S3 PUT price: $0.98
#Runtime: 9 hours 20 mins
println(path2cube);
dc = Zarr.zopen(path2cube)
@time for i = 1:n; dc["v"][r[i],c[i],:]; end


#=



path2cube = "http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr";
println(path2cube);
dc = Zarr.zopen(path2cube)
@time for i = 1:n; dc["v"][r[i],c[i],:]; end


path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_t_25730_xy_10_test.zarr";
println(path2cube)
dc = Zarr.zopen(path2cube);
@time for i = 1:n; dc["v"][r[i],c[i],:]; end

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_t_25730_xy_5_test.zarr";
println(path2cube)
dc = Zarr.zopen(path2cube);
@time for i = 1:n; dc["v"][r[i],c[i],:]; end

path2cube = "http://its-live-data.s3.amazonaws.com/test_datacubes/forAlex/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000_rechunked_t_25730_xy_1_test.zarr";
println(path2cube)
dc = Zarr.zopen(path2cube);
@time for i = 1:n; dc["v"][r[i],c[i],:]; end
=#