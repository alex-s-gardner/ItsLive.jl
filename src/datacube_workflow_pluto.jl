### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 214a6be6-dbe6-4281-97ef-ab5fd35b7e50
# import package manager and activate ITS_LIVE
begin
	import Pkg
	Pkg.activate(".")
	Pkg.activate("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl/Project.toml")
	include("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl")
	using ITS_LIVE, Zarr
end

# ╔═╡ 03d598ee-a11a-473f-9352-d9b813bddaa5
# load in ITS_LIVE datacube catalog as a Julia DataFrame
catalogdf = ITS_LIVE.catalog();

# ╔═╡ e6c9da34-c5fd-40a8-9a91-f406a64467c1
# find the DataFrame rows of the datacube that intersect a series of lat/lon points
lat = [69.1, 70.1, 90]

# ╔═╡ 049da617-6e50-40ee-80d4-36ce563e898a
lon = [-49.4, -48, 15]

# ╔═╡ 2723fb71-1687-4875-bb76-18ff8be7e50e
rows = ITS_LIVE.intersect.(lat,lon, Ref(catalogdf))

# ╔═╡ 00c7cebb-b314-4f37-8040-da6b405a9100
# find valid intersecting points
valid_intersect = .~ismissing.(rows)

# ╔═╡ 4c6970fd-a4ed-4bb3-a076-510a122bf361
begin
	# remove non-intersecting points 
	deleteat!(lat,findall(.~valid_intersect))
	deleteat!(lon,findall(.~valid_intersect))
	deleteat!(rows,findall(.~valid_intersect))
end

# ╔═╡ 9d875050-77f3-497a-bb73-4ac73e2dede9
# find all unique datacubes (mutiple points can intersect the same cube)
urows = unique(rows)

# ╔═╡ 7ffeabd1-4c09-423c-906b-6665b6197136
for row = urows
    # extract path to Zarr datacube file
    path2cube = catalogdf[row,"zarr_url"]

    # load Zarr group
    dc = Zarr.zopen(path2cube)

    # find closest point
   xind, yind = ITS_LIVE.nearestxy(lat, lon, dc)
end

# ╔═╡ c9b759d5-612a-4fe0-93dc-d617a2bda309


# ╔═╡ Cell order:
# ╠═214a6be6-dbe6-4281-97ef-ab5fd35b7e50
# ╠═03d598ee-a11a-473f-9352-d9b813bddaa5
# ╠═e6c9da34-c5fd-40a8-9a91-f406a64467c1
# ╠═049da617-6e50-40ee-80d4-36ce563e898a
# ╠═2723fb71-1687-4875-bb76-18ff8be7e50e
# ╠═00c7cebb-b314-4f37-8040-da6b405a9100
# ╠═4c6970fd-a4ed-4bb3-a076-510a122bf361
# ╠═9d875050-77f3-497a-bb73-4ac73e2dede9
# ╠═7ffeabd1-4c09-423c-906b-6665b6197136
# ╠═c9b759d5-612a-4fe0-93dc-d617a2bda309
