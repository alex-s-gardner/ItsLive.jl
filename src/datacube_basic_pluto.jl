### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 214a6be6-dbe6-4281-97ef-ab5fd35b7e50
begin
	md"""
	**Load packages, activate ITS_LIVE [~30s]:**
	"""
	import Pkg
	Pkg.activate(".")
	Pkg.activate("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl/Project.toml")
	include("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl")
	using ITS_LIVE, Zarr, Plots, Dates, DateFormats, PlutoUI
	Plots.PyPlotBackend()
end

# ╔═╡ a0982a9c-2436-4fed-bd18-6b2de4e4758c
md"""

##### ITS_LIVE: interactive glacier velocity time series notebook

A simple Pluto notebook showing how to use the basic functionality of the ITS\_LIVE.jl package. E.g. how to access and plot ITS_LIVE timeseries of glacier velocities for single lat/lon locations.

"""

# ╔═╡ 03d598ee-a11a-473f-9352-d9b813bddaa5
# load in ITS_LIVE datacube catalog as a Julia DataFrame
# [~5s]
catalogdf = ITS_LIVE.catalog();

# ╔═╡ f0a6d6d5-6557-4585-a176-e58494d85758
begin
	latin = @bind lat NumberField(-90:0.0001:90, default=60.0480)
	lonin = @bind lon NumberField(-180:0.0001:180, default=-140.5153)
	varin = @bind var Select(["v", "v_error"])
	
	md"""
	**Input latitude and longitude for point of interest [default = Malispina]:**
	
	Latitude: $(latin) [decimal degrees]
	
	Longitude: $(lonin) [decimal degrees]

	Variable: $(varin)
	"""
end

# ╔═╡ e02c67c1-4e18-4738-a4bd-60264a5e017d
# retrieve and plot data
begin
	# retrieve data from datacube sitting in the AWS cloud
	C = ITS_LIVE.getvar(lat,lon,["mid_date", var], catalogdf)
	# exclude missing data
    valid =  .!ismissing.(C[1,var]);
	# plot data
    plot(C[1,"mid_date"][valid], C[1,var][valid], seriestype = :scatter); 	
end

# ╔═╡ Cell order:
# ╟─a0982a9c-2436-4fed-bd18-6b2de4e4758c
# ╟─214a6be6-dbe6-4281-97ef-ab5fd35b7e50
# ╠═03d598ee-a11a-473f-9352-d9b813bddaa5
# ╟─f0a6d6d5-6557-4585-a176-e58494d85758
# ╠═e02c67c1-4e18-4738-a4bd-60264a5e017d
