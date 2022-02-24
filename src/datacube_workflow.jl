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

# filter data and plot
Plots.PyPlotBackend()
plot()
outlier = Vector{BitVector}()
dtmax = Vector{Vector{Union{Missing, Float64}}}()

for i = 1:lastindex(lat)

    outlie0, dtmax0, sensorgroups = ITS_LIVE.vxvyfilter(C[i,"vx"],C[i,"vy"],C[i,"date_dt"], C[i,"satellite_img1"])
    push!(outlier, outlie0)
    push!(dtmax, dtmax0)
    
    v = sqrt.(float.(C[i,"vx"]).^2 + float.(C[i,"vy"]).^2)

    if any(outlier[i])
        p = plot!(C[i,"mid_date"][outlier[i]], v[outlier[i]], seriestype = :scatter, mc = :gray)
    end

    valid =  .!ismissing.(v) .& .!outlier[i]
    p = plot!(C[i,"mid_date"][valid], v[valid], seriestype = :scatter)
    
    display(p)
end

# fit seasonal model with iterannual changes in amplitude and phase]
i = 4;
valid = (.!ismissing.(C[i,"vx"])) .& (.!outlier[i])
tx_fit, vx_fit, ampx_fit, phasex_fit, vx_fit_err, ampx_fit_err, fitx_count, fitx_outlier_frac, outlier[i][valid] = ITS_LIVE.lsqfit(C[i,"vx"][valid],C[i,"vx_error"][valid],C[i,"mid_date"][valid],C[i,"date_dt"][valid])

# interpolate lsqfit model to fine time resolution
t_i = DateTime.(DateFormats.YearDecimal.(2013:0.1:2022))
vx_i, vx_i_err = ITS_LIVE.lsqfit_interp(tx_fit, vx_fit, ampx_fit, phasex_fit, vx_fit_err, ampx_fit_err, t_i) # need to fix extraploaltion at some point

# plot data and fit
p = plot(C[i,"mid_date"][outlier[i]], C[i,"vx"][outlier[i]], seriestype = :scatter, mc = :gray)
p = plot!(C[i,"mid_date"][.!outlier[i]], C[i,"vx"][.!outlier[i]], seriestype = :scatter)
p = plot!(t_i, vx_i)