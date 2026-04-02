"""
    plot_disaggregate_pixel(vx, vy, t1, t2, vx_fit, vy_fit, valid_obs; kwargs...) -> Figure

Two-panel figure for a single ITS_LIVE pixel: observation segments coloured by
interval duration on the left, a zoomed detail panel on the right, and the
disaggregated speed signal ± 2σ overlaid on both panels.

# Arguments
- `vx`, `vy`: raw velocity components (m/yr); may contain `missing`
- `t1`, `t2`: observation interval start/end dates
- `vx_fit`, `vy_fit`: `DimStack`s returned by `ItsLive.disaggregate` (`:signal`, `:std`)
- `valid_obs`: `BitVector` — `true` = observation used in final fit

# Keyword Arguments
- `method`: disaggregation method (used for legend label); default `nothing` → `"fit"`
- `pt_name`: location name shown in the title (default `""`)
- `latlon`: `(lat, lon)` tuple shown in the title (default `(NaN, NaN)`)
- `detail_interval`: `(start_date, end_date)` for the zoom panel
  (default `(Date(2020,6,1), Date(2022,6,1))`)
"""
function plot_disaggregate_pixel(
    vx, vy, t1, t2,
    vx_fit, vy_fit,
    valid_obs;
    method = nothing, 
    pt_name::String = "",
    latlon = (NaN, NaN),
    detail_interval = (Date(2020, 6, 1), Date(2022, 6, 1)),
)

    method_name = isnothing(method) ? "fit" : nameof(typeof(method)) 

    output_start = first(dims(vx_fit, :Ti))
    output_end = last(dims(vx_fit, :Ti))

    not_missing = .!ismissing.(vx)
    v      = hypot.(vx[not_missing], vy[not_missing])
    pt1_pts = GI.Point.(yeardecimal.(t1[not_missing]), v)
    pt2_pts = GI.Point.(yeardecimal.(t2[not_missing]), v)

    dur     = Float32.((yeardecimal.(t2[not_missing]) .- yeardecimal.(t1[not_missing])) .* 365.25)
    dur_max = sort(dur)[end - max(0, div(length(dur), 20))]

    segs        = [(Point2f(p1.geom...), Point2f(p2.geom...)) for (p1, p2) in zip(pt1_pts, pt2_pts)]
    seg_colors  = repeat(dur; inner=2)
    seg_cmap    = reverse(sequential_palette(300, 256; c=0.88, s=0.6))

    col_flt2_sig = RGBAf(0.95, 0.38, 0.00, 1.00)
    col_rem      = RGBAf(0.80, 0.80, 0.80, 0.90)

    segs_removed = segs[.!valid_obs[not_missing]]

    zoom_xlims = (yeardecimal(detail_interval[1]), yeardecimal(detail_interval[2]))
    obs_mask   = yeardecimal.(t1[not_missing]) .>= zoom_xlims[1] .&&
                 yeardecimal.(t2[not_missing]) .<= zoom_xlims[2]

    decyear_out = yeardecimal.(dims(vx_fit, Ti).val)
    out_mask    = decyear_out .>= zoom_xlims[1] .&& decyear_out .<= zoom_xlims[2]

    v2     = hypot.(vx_fit[:signal].data, vy_fit[:signal].data)
    v2_err = sqrt.((vx_fit[:std].data .* abs.(vx_fit[:signal].data)).^2 .+
                   (vy_fit[:std].data .* abs.(vy_fit[:signal].data)).^2) ./ max.(v2, 1e-6)

    ax1_mask = yeardecimal.(t1[not_missing]) .>= minimum(decyear_out) .&&
               yeardecimal.(t2[not_missing]) .<= maximum(decyear_out)
    ax1_ymax = maximum(vcat(Float64.(v[ax1_mask]), v2 .+ 2 .* v2_err))
    ax1_ypad = 0.05 * ax1_ymax

    zoom_yvals = vcat(collect(extrema(v[obs_mask])), collect(extrema(v2[out_mask])))
    zoom_ymin  = minimum(zoom_yvals)
    zoom_ymax  = maximum(zoom_yvals)
    zoom_ypad  = 0.1 * (zoom_ymax - zoom_ymin)
    zoom_ylims = (max(0.0, zoom_ymin - zoom_ypad), zoom_ymax + zoom_ypad)

    yr_start = year(output_start)
    yr_end   = ceil(Int, Float64(year(output_end))) + 1

    fig = Figure(size=(1380, 480), fontsize=12,
        fonts=(; regular="Helvetica", bold="Helvetica Bold"))

    lat_str = isnan(latlon[1]) ? "" :
        "\n$(round(latlon[1], digits=3))°N, $(round(abs(latlon[2]), digits=3))°W"

    ax1 = Axis(fig[1, 1];
        ylabel    = "Speed (m yr⁻¹)",
        title     = "Single-pixel disaggregation — $pt_name$lat_str",
        titlesize = 12,
        xticks    = Float64(yr_start):2.0:Float64(yr_end),
    )

    ax2 = Axis(fig[1, 2];
        title              = "$(Dates.format(detail_interval[1], "yyyy")) – " *
                             "$(Dates.format(detail_interval[2], "yyyy")) detail",
        titlesize          = 12,
        yticklabelsvisible = false,
        yticksvisible      = false,
        xticks             = floor(zoom_xlims[1]):1.0:ceil(zoom_xlims[2]),
        xminorticks        = floor(zoom_xlims[1]):(1/12):ceil(zoom_xlims[2]),
        xminorticksvisible = true,
        xminorticksize     = 3,
    )

    for ax in (ax1, ax2)
        for yr in yr_start:yr_end
            isodd(yr) && vspan!(ax, Float64(yr), Float64(yr+1); color=RGBA(0,0,0,0.025))
        end
        linesegments!(ax, segs;
            color=seg_colors, colormap=seg_cmap,
            colorrange=(0f0, Float32(dur_max)), linewidth=0.8)
        linesegments!(ax, segs_removed; color=col_rem, linewidth=0.8, label="outlier")
        band!(ax, decyear_out, v2 .- 2 .* v2_err, v2 .+ 2 .* v2_err;
            color=(col_flt2_sig, 0.20))
        lines!(ax, decyear_out, v2;
            color=col_flt2_sig, linewidth=2.5,
            label="$(method_name) ± 2σ")
    end

    xlims!(ax1, minimum(decyear_out), maximum(decyear_out))
    ylims!(ax1, 0.0, ax1_ymax + ax1_ypad)
    xlims!(ax2, zoom_xlims...)
    ylims!(ax2, zoom_ylims...)

    lines!(ax1,
        [zoom_xlims[1], zoom_xlims[2], zoom_xlims[2], zoom_xlims[1], zoom_xlims[1]],
        [zoom_ylims[1], zoom_ylims[1], zoom_ylims[2], zoom_ylims[2], zoom_ylims[1]];
        color=RGBA(0,0,0,0.55), linewidth=1.4, linestyle=:dot)

    axislegend(ax1; position=:lt, framevisible=true,
        framecolor=RGBA(0,0,0,0.15), labelsize=10, rowgap=3)

    Label(fig[2, 1:2]; text="Year", fontsize=12, tellwidth=false)

    Colorbar(fig[1, 3];
        colormap=seg_cmap, colorrange=(0f0, Float32(dur_max)),
        label="Observation interval (days)",
        labelsize=11, ticksize=4, width=14)

    colsize!(fig.layout, 1, Relative(0.60))
    colsize!(fig.layout, 2, Relative(0.30))
    rowsize!(fig.layout, 2, Fixed(22))
    colgap!(fig.layout, 1, 30)
    colgap!(fig.layout, 2, 8)
    rowgap!(fig.layout, 1, 4)

    return fig
end