# ITS_LIVE.jl
![ITS_LIVE](https://its-live-data.s3.amazonaws.com/documentation/ITS_LIVE_logo_transparent_wht.png)

# A Julia Package 
**This repository is the beginnings of a Julia package for working with NASA [ITS_LIVE](https://its-live.jpl.nasa.gov/) data, it is in its infancy and will be developed over time**

## Function List 
**`catalog`** returns a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes

**`intersect`** returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longituded

**`nearestxy`** function returns the x/y indicies into a ZarrGroup for the points nearest the provided lat, lon locations
