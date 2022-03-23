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
# import package manager, load packages, activate ITS_LIVE, and set plotting backend
# [~30s]
begin
	import Pkg
	Pkg.activate(".")
	Pkg.activate("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl/Project.toml")
	include("/Users/gardnera/Documents/GitHub/ITS_LIVE.jl")
	using ITS_LIVE, Zarr, Plots, Dates, DateFormats, PlutoUI
	Plots.PyPlotBackend()
end

# ╔═╡ 03d598ee-a11a-473f-9352-d9b813bddaa5
# load in ITS_LIVE datacube catalog as a Julia DataFrame
# [~5s]
catalogdf = ITS_LIVE.catalog();

# ╔═╡ ce4eb188-d8e1-4214-b0e2-5b3301e4ab38
@bind lat NumberField(-90:0.0001:90, default=60.048110121383285)

# ╔═╡ f0a6d6d5-6557-4585-a176-e58494d85758
@bind lon NumberField(-180:0.0001:180, default=-140.5153002286824)

# ╔═╡ e6c9da34-c5fd-40a8-9a91-f406a64467c1
# define points of interest malispina example
#begin
#	lat = [59.92518849057908, 60.00020529502237, 60.048010121383285];
#	lon = [-140.62047084667643, -140.5638405139104, -140.5153002286824];
#end

# ╔═╡ 2723fb71-1687-4875-bb76-18ff8be7e50e
# retrieve data columns from Zarr as a named matrix of vectors
# [~17s]
begin
	varnames = ["mid_date", "v"];
	C = ITS_LIVE.getvar(lat,lon,varnames, catalogdf);
end

# ╔═╡ e02c67c1-4e18-4738-a4bd-60264a5e017d
# plot data
begin
	plot()
    valid =  .!ismissing.(C[i,"vx"]);
    plot!(C[i,"mid_date"][valid], C[i,"v"][valid], seriestype = :scatter);
    	
end
