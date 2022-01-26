trans = Proj4.Transformation("EPSG:4326", "EPSG:" * dc.attrs["projection"])
x, y = trans([-40,80])x