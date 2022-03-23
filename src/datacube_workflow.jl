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

# noisy data test
# lat = [60.6273]
# lon = [-140.1370]

# slow moving test
# lat = [59.823]
# lon= [-138.179]

# malispina example
lat = [59.92518849057908, 60.00020529502237, 60.048010121383285, 60.08331214700366, 60.12523330242785, 60.168625375586224, 60.21348836647878, 60.25614498077006]
lon = [-140.62047084667643, -140.5638405139104, -140.5153002286824, -140.46749540232148, -140.43881250650492, -140.43219337670112, -140.38733038580855, -140.33143551190963]

# jakobshavn glacier example
#lat = [69.1302, 69.1155, 69.1020, 69.0994, 69.1034, 69.1074, 69.1101, 69.1148]
#lon = [-49.5229, -49.4833, -49.4430, -49.3987, -49.3463, -49.3080, -49.2691, -49.2267]

# example of datacube contents
# dc = Zarr.zopen(catalogdf[1,"zarr_url"]).arrays

# retrieve data columns from Zarr as a named matrix of vectors
varnames = ["mid_date", "date_dt", "vx", "vx_error", "vy", "vy_error","satellite_img1"]
@time C = ITS_LIVE.getvar(lat,lon,varnames, catalogdf)

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
i = 1;|
valid = (.!ismissing.(C[i,"vx"])) .& (.!outlier[i])
model = "sinusoidal"
tx_fit, vx_fit, vx_amp, vx_phase, vx_fit_err, vx_amp_err, vx_fit_count, vx_fit_outlier_frac, outlier[i][valid] = 
    ITS_LIVE.lsqfit(C[i,"vx"][valid],C[i,"vx_error"][valid],C[i,"mid_date"][valid],C[i,"date_dt"][valid]; model = model);

ty_fit, vy_fit, vy_amp, vy_phase, vy_fit_err, vy_amp_err, vy_fit_count, vy_fit_outlier_frac, outlier[i][valid] = 
    ITS_LIVE.lsqfit(C[i,"vy"][valid],C[i,"vy_error"][valid],C[i,"mid_date"][valid],C[i,"date_dt"][valid]; model = model);

# solve for velocity in in 2017.5 and velocity trend
vx0, dvx_dt, vx0_err = ITS_LIVE.wlinearfit(tx_fit, vx_fit, vx_fit_err, DateTime(2017,7,1))
vy0, dvy_dt, vy0_err = ITS_LIVE.wlinearfit(ty_fit, vy_fit, vy_fit_err, DateTime(2017,7,1))

# compute velocity magnitude metrics from component values
v_fit, v_fit_err, v_fit_count, v_fit_outlier_frac  = ITS_LIVE.annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, vy_fit_err, vx_fit_count, vy_fit_count, vx_fit_outlier_frac, vy_fit_outlier_frac)
v0, v0_err, dv0_dt, v0_amp, v0_amp_err, v0_phase = ITS_LIVE.climatology_magnitude(vx0, vy0, vx0_err, vy0_err, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err, vx_phase, vy_phase)

# interpoalte model in time 
t_i = DateTime.(DateFormats.YearDecimal.(2013:0.1:2022))
vx_i, vx_i_err = ITS_LIVE.lsqfit_interp(tx_fit, vx_fit, vx_amp, vx_phase, vx_fit_err, vx_amp_err, t_i; interp_method = "BSpline"); # need to fix extraploaltion at some point
vy_i, vy_i_err = ITS_LIVE.lsqfit_interp(ty_fit, vy_fit, vy_amp, vy_phase, vy_fit_err, vy_amp_err, t_i; interp_method = "BSpline"); # need to fix extraploaltion at some point

# plot data and fit
p = plot(C[i,"mid_date"][outlier[i]], C[i,"vx"][outlier[i]], seriestype = :scatter, mc = :gray)
p = plot!(C[i,"mid_date"][.!outlier[i]], C[i,"vx"][.!outlier[i]], seriestype = :scatter)
p = plot!(t_i, vx_i)



# write data to a .csv
#=
using DelimitedFiles

open("/Users/gardnera/Downloads/MalPt.csv"; write=true) do f
    write(f, "rows = " * "mid_date - date_dt - vx - vx_error - vy - vy_error - satellite_img1" * ", \n")
    writedlm(f,  C[3:end], ',')
end
=#