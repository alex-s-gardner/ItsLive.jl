# ITS_LIVE.jl
![ITS_LIVE](https://its-live-data.s3.amazonaws.com/documentation/ITS_LIVE_logo_transparent_wht.png)

# A Julia Package 
**This repository is the beginnings of a Julia package for working with NASA [ITS_LIVE](https://its-live.jpl.nasa.gov/) data, it is in its infancy and will be developed over time**

## Working Example
**`datacube_workflow.jl`** example script showing how to work with the `ITS_LIVE.jl` package. `datacube_workflow.jl` should run smoothly after locally cloning `ITS_LIVE.jl`.

**`datacube_workflow_pluto.jl`** [NOT FINISHED] Pluto implementation of `datacube_workflow.jl`.

## Function List 
**`catalogdf = ITS_LIVE.catalog()`** returns a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes

**`rownumdf = ITS_LIVE.intersect(latitude, longitude, catalogdf)`** returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longitude

**`x, y = ITS_LIVE.nearestxy(latitude, longitude, DataArray)`** returns the x/y indices into a ZarrGroup for the points nearest the provided lat, lon locations

**`M = ITS_LIVE.getvar(latitude, longitude, varnames, catalogdf)`** returns a named m x n matrix of vectors with m = length(lat) rows and n = length(varnames)+2(for lat and lon) columns for the points nearest the lat/lon location from ITS_LIVE Zarr datacubes

**`dtmax = dtfilter(x,dt)`** returns the maximum `dt` for which the distribution of x shows no statistical difference from the distribution of x in the minimum `dt` bin

This filter is needed to identify longer dts that exhibit "skipping" or "locking" behavior in feature tracking estimates of surface flow. This happens when the surface texture provides a lesser match than to stationary features, due to long time separation between repeat images, such as ice falls and curved medial moraines.

**`outlier, dtmax, sensorgroups = vxvyfilter(x,dt)`** applies `dtfilter()` to `vx` and `vy` and returns `outlier` BitVector, the maximum `dt` (`dtmax`) for which the distribution of vx *or* vy shows no statistical difference and the sensor groupings used when filtering the data (`sensorgroups`)

**`t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, fit_count, fit_outlier_frac, outlier = lsqfit(v,v_err,mid_date,date_dt,mad_thresh)`** error wighted model fit to discrete interval data. The current model is an iterative fit to a sinusoidal function with unique mean (`v_fit`), phase (`phase_fit`) and amplitude (`amp_fit`) for each year of data centered at time `t_fit`. Errors are provided for annual means (`v_fit_err`) and amplitude (`amp_fit_err`). `fit_count` gives the number of velocity observations used for each year,  `fit_outlier_frac` provides the fraction of data excluded from the model fit and a `outlier` BitVector 

**`v_i, v_i_err = lsqfit_interp(t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, t_i)`** evaluates the outputs of `lsqfit` at times `t_i`, outputting velocity (`v_i`) and velocity error (`v_i_err`) at times `t_i`.

**`v0 = running_mean(v, w)`** calculates the running mean of `v` with window size of `w`

**`offset, slope, error = wliearfit(t, v, v_err, datetime0)`** returns the `offset`, `slope`, and  error (`error`) for a weighted linear fit to `v` with a y-intercept of `datetime0`

 **`v, v,_err, dv_dt, v_amp, v_amp_err, v_phase = climatology_magnitude(vx0, vy0, vx0_err, vy0_err, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err, vx_phase, vy_phase)`** returns the mean (`v`), standard error (`v_err`), trend (`dv_dt`), seasonal amplitude (`v_amp`), error in seasonal amplitude (`v_amp_err`), and seasonal phase (`v_phase`) from component values projected on the unit flow vector defined by vx0 and vy0

  **`v_fit, v_fit_err, v_fit_count, v_fit_outlier_frac  = annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, vy_fit_err, vx_fit_count, vy_fit_count, vx_fit_outlier_frac, vy_fit_outlier_frac)`** returns the annual mean (`v_fit`), error (`v_fit_err`), count(`v_fit_count`), and outlier fraction (`v_fit_outlier_frac`) from component values projected on the unit flow  vector defined by vx0 and vy0

