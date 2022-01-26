# ITS_LIVE.jl
![ITS_LIVE](https://its-live-data.s3.amazonaws.com/documentation/ITS_LIVE_logo_transparent_wht.png)

# A Package Julia
**This repository is the beginnings of a Julia package for working with NASA [ITS_LIVE](https://its-live.jpl.nasa.gov/) data, it is in it's infancy and will be developed over time**

## Function List 
**`datacube_catalog`** returns a DataFrame of the catalog for all of the ITS_LIVE zarr datacubes
**`datacube_catalog`** returns the rownumber of the the DataFrame catalog of the ITS_LIVE zarr datacubes that intersects the provided latitude and longituded
**`datacube_nearestxy`** function returns the x/y indicies into a ZarrGroup for the points nearest the provided lat, lon locations
