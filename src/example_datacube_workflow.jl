"""
    this is an example workflow for working with ITS_LIVE.jl and the zarr datacubes

# Author
Alex S. Gardner, JPL, Caltech.
"""
# Revise.jl allows you to modify code and use the changes without restarting Julia
using Plots, Dates, DateFormats, ItsLive, HDF5

# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = ItsLive.catalog()

# find the DataFrame rows of the datacube that intersect a series of lat/lon points

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
@time C = ItsLive.getvar(lat,lon,varnames, catalogdf)

# plot a variable for all points
ItsLive.plotvar(C,"vx")

# for a single point do some filtering
i = 1;

# plot by sensor
p = ItsLive.plotbysensor(C[i,:],"vx"); plot(p)

# filter long dt data that exhibits skipping behaviour
outlier, dtmax, sensorgroups = ItsLive.vxvyfilter(C[i,"vx"],C[i,"vy"],C[i,"date_dt"]; sensor = C[i,"satellite_img1"])

# fit seasonal model with iterannual changes in amplitude and phase]
valid = (.!ismissing.(C[i,"vx"])) .& (.!outlier)
model = "sinusoidal"
tx_fit, vx_fit, vx_amp, vx_phase, vx_fit_err, vx_amp_err, vx_fit_count, vx_fit_outlier_frac, outlier[valid] = 
    ItsLive.lsqfit_annual(C[i,"vx"][valid],C[i,"vx_error"][valid],C[i,"mid_date"][valid],C[i,"date_dt"][valid]; model = model);

# interpoalte model in time 
t_i = DateTime.(DateFormats.YearDecimal.(2013:0.1:2022))
vx_i, vx_i_err = ItsLive.lsqfit_interp(tx_fit, vx_fit, vx_amp, vx_phase, vx_fit_err, vx_amp_err, t_i; interp_method = "BSpline"); # need to fix extraploaltion at some point

# plot data and fit
plot()
p = plot!(C[i,"mid_date"][.!outlier], C[i,"vx"][.!outlier], seriestype = :scatter)
p = plot!(t_i, vx_i)
