# ITS_LIVE.jl
![ITS_LIVE](https://its-live-data.s3.amazonaws.com/documentation/ITS_LIVE_logo_transparent_wht.png)

# A Julia Package 
**This repository is the beginnings of a Julia package for working with NASA [ITS_LIVE](https://its-live.jpl.nasa.gov/) data, it is in its infancy and will be developed over time**

## Function List 
**`catalogdf = ITS_LIVE.catalog()`** returns a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes

**`rownumdf = ITS_LIVE.intersect(latitude, longitude, catalogdf)`** returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longitude

**`x, y = ITS_LIVE.nearestxy(latitude, longitude, DataArray)`** returns the x/y indices into a ZarrGroup for the points nearest the provided lat, lon locations

**`M = ITS_LIVE.getvar(latitude, longitude, varnames, catalogdf)`** returns a named m x n matrix of vectors with m = length(lat) rows and n = length(varnames)+2(for lat and lon) columns for the points nearest the lat/lon location from ITS_LIVE Zarr datacubes

**`dtmax = dtfilter(x,dt)`** returns the maximum `dt` for which the distribution of x shows no statistical difference from the distribution of x in the minimum `dt` bin

This filter is needed to identify longer dts that exhibit "skipping" or "locking" behavior in feature tracking estimates of surface flow. This happens when the surface texture provides a lesser match than to stationary features, due to long time separation between repeat images, such as ice falls and curved medial moraines.

**`vx, vy, dtmax, outlierfrac, sensorgoup = vxvyfilter(x,dt)`** applies `dtfilter()` to `vx` and `vy` and returns the maximum `dt` (`dtmax`) for which the distribution of vx *or* vy shows no statistical difference and the fraction of data that was removed (`outlierfrac`)

