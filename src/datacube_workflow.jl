"""
    this is an example workflow for working with ITS_LIVE.jl and the zarr datacubes

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 25, 2022
"""

# Revise.jl allows you to modify code and use the changes without restarting Julia
using Zarr, Revise, Plots, Dates, DateFormats

import ITS_LIVE

# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = ITS_LIVE.catalog()

# find the DataFrame rows of the datacube that intersect a series of lat/lon points
lat = [69.1302, 69.1155, 69.1020, 69.0994, 69.1034, 69.1074, 69.1101, 69.1148]
lon = [-49.5229, -49.4833, -49.4430, -49.3987, -49.3463, -49.3080, -49.2691, -49.2267]

# example of datacube contents
dc = Zarr.zopen(catalogdf[1,"zarr_url"]).arrays

# retrieve data columns from Zarr as a named matrix of vectors
varnames = ["mid_date", "date_dt", "vx", "vx_error", "vy", "vy_error","satellite_img1"]
@time C = ITS_LIVE.getvar(lat,lon,varnames ,catalogdf)

# make same plot but filter data using vxvyfilter
Plots.PyPlotBackend()
plot()
for i = 1:lastindex(lat)
    C[i,"vx"], C[i,"vy"], dtmax, outlierfrac, sensorgroups = ITS_LIVE.vxvyfilter(C[i,"vx"],C[i,"vy"],C[i,"date_dt"], C[i,"satellite_img1"])
    v = sqrt.(float.(C[i,"vx"]).^2 + float.(C[i,"vy"]).^2)
    p = plot!(C[i,"mid_date"], v, seriestype = :scatter)
    display(p)
end

# apply a least squares model fit 
i = 1;

ti = DateTime.(DateFormats.YearDecimal.(2014:0.01:2021))

@time for k = 1:100
    ampx,phasex,amp_errx,t_intx,v_intx,v_int_errx,N_intx,outlier_fracx = ITS_LIVE.lsqfit(C[i,"vx"],C[i,"vx_error"],C[i,"mid_date"],C[i,"date_dt"])

    vix, v_int_ix = lsqfit_interp(t_intx,v_intx,ampx,phasex,ti,v_int_errx,amp_errx) # need to fix extraploaltion at some point

    p = plot(C[i,"mid_date"], C[i,"vx"], seriestype = :scatter)
    p = plot!(ti, vix)


    ampy,phasey,amp_erry,t_inty,v_inty,v_int_erry,N_inty,outlier_fracy = ITS_LIVE.lsqfit(C[i,"vy"],C[i,"vy_error"],C[i,"mid_date"],C[i,"date_dt"])
end

v_int = sqrt.(v_intx.^2 + v_inty.^2)
plot!(t_intx,v_int,seriestype = :scatter, color = :red)