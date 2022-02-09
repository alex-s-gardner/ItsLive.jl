"""
    this is an example workflow for working with ITS_LIVE zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 25, 2022
"""

# Revise.jl allows you to modify code and use the changes without restarting Julia
using Zarr, Revise, Plots

import ITS_LIVE

# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = ITS_LIVE.catalog()

# find the DataFrame rows of the datacube that intersect a series of lat/lon points
lat = [69.1, 70.1, 90]
lon = [-49.4, -48, 15]

# example of datacube contents
dc = Zarr.zopen(catalogdf[1,"zarr_url"]).arrays

# retrieve data columns from Zarr as a matrix of vectors
varnames = ["mid_date", "date_dt", "vx", "vx_error", "vy", "vy_error"]
@time M = ITS_LIVE.getvar(lat,lon,varnames ,catalogdf)

plot(M[1,"mid_date"], M[1,"vx"], seriestype = :scatter)