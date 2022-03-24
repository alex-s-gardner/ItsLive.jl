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
# import package manager, load packages, activate ItsLive, and set plotting backend
# [~30s]
begin
	import Pkg
	Pkg.activate(".")
	Pkg.activate("/Users/gardnera/Documents/GitHub/ItsLive.jl/Project.toml")
	include("/Users/gardnera/Documents/GitHub/ItsLive.jl")
	using ItsLive, Zarr, Plots, Dates, DateFormats, PlutoUI
	Plots.PyPlotBackend()
end

# ╔═╡ 03d598ee-a11a-473f-9352-d9b813bddaa5
# load in ITS_LIVE datacube catalog as a Julia DataFrame
# [~5s]
catalogdf = ItsLive.catalog();

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
	varnames = ["mid_date", "date_dt", "vx", "vx_error", 
		"vy", "vy_error","satellite_img1"];
	C = ItsLive.getvar(lat,lon,varnames, catalogdf);
end

# ╔═╡ a75a94b9-c143-4330-b644-4315782aaf22
# initialize arrays
begin
	outlier = Vector{BitVector}();
	dtmax = Vector{Vector{Union{Missing, Float64}}}();
end

# ╔═╡ fe195d74-6f77-4559-8e07-7a7db4f90800
# filter data based on time seperation between image-pairs
begin
	i = 1
	outlie0, dtmax0, sensorgroups = 
		ItsLive.vxvyfilter(C[i,"vx"],C[i,"vy"],C[i,"date_dt"], 
			C[i,"satellite_img1"]);
    push!(outlier, outlie0);
    push!(dtmax, dtmax0);
end

# ╔═╡ e02c67c1-4e18-4738-a4bd-60264a5e017d
begin
	plot()

	# calculate velocity magnitude
    vv = sqrt.(float.(C[i,"vx"]).^2 + float.(C[i,"vy"]).^2);

    if any(outlier[i])
        plot!(C[i,"mid_date"][outlier[i]], 
			vv[outlier[i]], seriestype = :scatter, mc = :gray);
    end

    valid =  .!ismissing.(vv) .& .!outlier[i];
    plot!(C[i,"mid_date"][valid], vv[valid], seriestype = :scatter);
    	
end

# ╔═╡ 04588ca9-9cd6-4dcb-8459-3009f4e2e818
# fit seasonal model with iterannual changes in amplitude and phase]
begin
	# valid = (.!ismissing.(C[i,"vx"])) .& (.!outlier[i])
	model = "sinusoidal";
	tx_fit, vx_fit, vx_amp, vx_phase, vx_fit_err, vx_amp_err, vx_fit_count,
	vx_fit_outlier_frac, outlier[i][valid] = 
		ItsLive.lsqfit(C[i,"vx"][valid],C[i,"vx_error"][valid],
			C[i,"mid_date"][valid],C[i,"date_dt"][valid]; model = model);
	
	ty_fit, vy_fit, vy_amp, vy_phase, vy_fit_err, vy_amp_err, vy_fit_count,
	vy_fit_outlier_frac, outlier[i][valid] = 
		ItsLive.lsqfit(C[i,"vy"][valid],C[i,"vy_error"][valid],
			C[i,"mid_date"][valid],C[i,"date_dt"][valid]; model = model);
end

# ╔═╡ c830c54b-f7e4-4ed0-acdf-9787dc08fd40
# solve for velocity in in 2017.5 and velocity trend
begin
	vx0, dvx_dt, vx0_se = ItsLive.wlinearfit(tx_fit, vx_fit, vx_fit_err,
		DateTime(2017,7,2));
	vy0, dvy_dt, vy0_se = ItsLive.wlinearfit(ty_fit, vy_fit, vy_fit_err,
		DateTime(2017,7,2));
end

# ╔═╡ 74e3bb68-df56-4fcb-a588-6fadf8afad52
# compute velocity magnitude metrics from component values
begin
	v_fit, v_fit_err, v_fit_count, v_fit_outlier_frac  = 
		ItsLive.annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, 
			vy_fit_err, vx_fit_count, vy_fit_count, vx_fit_outlier_frac, 
			vy_fit_outlier_frac);
	v, v_se, dv_dt, v_amp, v_amp_err, v_phase = ItsLive.climatology_magnitude(vx0,
		vy0, vx0_se, vy0_se, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err,
		vx_phase, vy_phase);
end

# ╔═╡ b73c7893-1601-4860-8cdd-74634f047f1b
# interpoalte model in time 
begin
	t_i = DateTime.(DateFormats.YearDecimal.(2013:0.1:2022))
	vx_i, vx_i_err = ItsLive.lsqfit_interp(tx_fit, vx_fit, vx_amp, vx_phase,
		vx_fit_err, vx_amp_err, t_i; interp_method = "BSpline"); 
	
	vy_i, vy_i_err = ItsLive.lsqfit_interp(ty_fit, vy_fit, vy_amp, vy_phase,
		vy_fit_err, vy_amp_err, t_i; interp_method = "BSpline"); 
end


# ╔═╡ 9d875050-77f3-497a-bb73-4ac73e2dede9
begin
	# plot data and fit
	plot(C[i,"mid_date"][outlier[i]], C[i,"vx"][outlier[i]], seriestype = :scatter, mc = :gray)
	plot!(C[i,"mid_date"][.!outlier[i]], C[i,"vx"][.!outlier[i]], seriestype = :scatter)
	plot!(t_i, vx_i)
end

# ╔═╡ c9b759d5-612a-4fe0-93dc-d617a2bda309


# ╔═╡ Cell order:
# ╠═214a6be6-dbe6-4281-97ef-ab5fd35b7e50
# ╠═03d598ee-a11a-473f-9352-d9b813bddaa5
# ╠═ce4eb188-d8e1-4214-b0e2-5b3301e4ab38
# ╠═f0a6d6d5-6557-4585-a176-e58494d85758
# ╠═e6c9da34-c5fd-40a8-9a91-f406a64467c1
# ╠═2723fb71-1687-4875-bb76-18ff8be7e50e
# ╠═a75a94b9-c143-4330-b644-4315782aaf22
# ╠═fe195d74-6f77-4559-8e07-7a7db4f90800
# ╠═e02c67c1-4e18-4738-a4bd-60264a5e017d
# ╠═04588ca9-9cd6-4dcb-8459-3009f4e2e818
# ╠═c830c54b-f7e4-4ed0-acdf-9787dc08fd40
# ╠═74e3bb68-df56-4fcb-a588-6fadf8afad52
# ╠═b73c7893-1601-4860-8cdd-74634f047f1b
# ╠═9d875050-77f3-497a-bb73-4ac73e2dede9
# ╠═c9b759d5-612a-4fe0-93dc-d617a2bda309
