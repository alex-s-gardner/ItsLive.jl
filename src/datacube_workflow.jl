"""
    this is an example workflow for working with ITS_LIVE zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 25, 2022
"""

# Revise.jl allows you to modify code and use the changes without restarting Julia
using Zarr, Revise

import ITS_LIVE

# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = ITS_LIVE.catalog()

# find the DataFrame rows of the datacube that intersect a series of lat/lon points
lat = [69.1, 70.1, 90]
lon = [-49.4, -48, 15]
rows = ITS_LIVE.intersect.(lat,lon, Ref(catalogdf))

# find valid intersecting points
valid_intersect = .~ismissing.(rows)

# remove non-intersecting points 
deleteat!(lat,findall(.~valid_intersect))
deleteat!(lon,findall(.~valid_intersect))
deleteat!(rows,findall(.~valid_intersect))

# find all unique datacubes (mutiple points can intersect the same cube)
urows = unique(rows)

yind = zeros(size(lat))
xind = zeros(size(lat))
for row = urows
    # extract path to Zarr datacube file
    path2cube = catalogdf[row,"zarr_url"]

    # load Zarr group
    dc = Zarr.zopen(path2cube)

    # find closest point
     a,b = ITS_LIVE.nearestxy(lat[row.==rows], lon[row.==rows], dc)
     xind[row.==rows] .= a
     yind[row.==rows] .= b

     # extract timeseries from datacube
     col1, col2 = [foo[i...,:] for i in ((1,2),(3,4))];

     mask = falses(size(dc["vx"]))
for 

     mask[[CartesianIndex((1,2)),CartesianIndex((1,2))],:] .= true
     foo = dc["vx"][[mask]


a = CartesianIndices((1,[1,4],1:size(dc["v"])[3]))

@time foo = 


mask = falses(4,5,1)
mask[3,2:4,1] .= true

mask


CartesianIndex(([1,3],[2,5]))
