# ITS_LIVE.jl
![ITS_LIVE](https://its-live-data.s3.amazonaws.com/documentation/ITS_LIVE_logo_transparent_wht.png)

# A Julia Package 
**This repository is the beginnings of a Julia package for working with NASA [ITS_LIVE](https://its-live.jpl.nasa.gov/) data, it is in its infancy and will be developed over time**

## Working Example
**`datacube_workflow.jl`** example script showing how to work with the `ITS_LIVE.jl` package. `datacube_workflow.jl` should run smoothly after locally cloning `ITS_LIVE.jl`.

## Function List 
**`catalogdf = ITS_LIVE.catalog()`** returns a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes

**`rownumdf = ITS_LIVE.intersect(latitude, longitude, catalogdf)`** returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longitude

**`x, y = ITS_LIVE.nearestxy(latitude, longitude, DataArray)`** returns the x/y indices into a ZarrGroup for the points nearest the provided lat, lon locations

**`M = ITS_LIVE.getvar(latitude, longitude, varnames, catalogdf)`** returns a named m x n matrix of vectors with m = length(lat) rows and n = length(varnames)+2(for lat and lon) columns for the points nearest the lat/lon location from ITS_LIVE Zarr datacubes

**`dtmax = dtfilter(x,dt)`** returns the maximum `dt` for which the distribution of x shows no statistical difference from the distribution of x in the minimum `dt` bin

This filter is needed to identify longer dts that exhibit "skipping" or "locking" behavior in feature tracking estimates of surface flow. This happens when the surface texture provides a lesser match than to stationary features, due to long time separation between repeat images, such as ice falls and curved medial moraines.

**`vx, vy, dtmax, outlierfrac, sensorgoup = vxvyfilter(x,dt)`** applies `dtfilter()` to `vx` and `vy` and returns the maximum `dt` (`dtmax`) for which the distribution of vx *or* vy shows no statistical difference and the fraction of data that was removed (`outlierfrac`)

**`t_fit, v_fit, amp_fit, phase_fit, amp_fit_err, v_fit_err, fit_count, fit_outlier_frac = lsqfit(v,v_err,mid_date,date_dt,mad_thresh)`** error wighted model fit to discrete interval data. The current model is an iterative fit to a sinusoidal function with unique mean (`v_fit`), phase (`phase_fit`) and amplitude (`amp_fit`) for each year of data centered at time `t_fit`. Errors are provided for amplitude (`amp_fit_err`) and annual means (`v_fit_err`). `fit_count` gives the number of velocity observations used for each year and  `fit_outlier_frac` provides the fraction of data excluded from the model fit. 

**`v_i, v_i_err = lsqfit_interp(t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, t_i)`** evaluates the outputs of `lsqfit` at times `t_i`, outputting velocity (`v_i`) and velocity error (`v_i_err`) at times `t_i`.


