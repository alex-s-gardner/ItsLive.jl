"""
    this is an example workflow for working with ITS_LIVE zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 25, 2022
"""

using AWS, Zarr, ArchGDAL, DataFrames

# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = datacube_catalog()

# find the DataFrame rows of the datacube that intersect a series of lat/lon points
lat = [69.1, 70.1, 90]
lon = [-49.4, -48, 15]
rows = datacube_intersect.(lat,lon, Ref(catalogdf))

# find valid intersecting points
valid_intersect = .~ismissing.(rows)

# remove non-intersecting points 
lat = lat[valid_intersect]
lon = lon[valid_intersect]
rows = rows[valid_intersect]

# find all unique datacubes (mutiple points can intersect the same cube)
urows = unique(rows)

for row = urows
    # extract path to Zarr datacube file
    path2cube = catalogdf[row,"zarr_url"]

    # load Zarr group
    dc = zopen(path2cube)

    # find closest point
   xind, yind = datacube_nearestxy(lat, lon, dc)
end
