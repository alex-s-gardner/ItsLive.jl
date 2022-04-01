"""
    p = plotvar(C::Named Matrix{Any}, varname::String)

    plot its_live data (`C`) variable (`varname`) for multiple points

# Example no inputs
```julia
julia> p = plotvar(C, varname)
```

# Arguments
   - `C::Named Matrix{Any}`: ITS_LIVE named matrix
   - `varname::String`: variable name to plot
# Keyword Arguments
    - `dtmax::Number`: maximum time seperation between image pairs [days] to plot

# Author
Alex S. Gardner, JPL, Caltech.
"""
function plotvar(C, varname::String; dtmax::Number = Inf)
    npts = size(C,1)
    Plots.PlotlyBackend()
    p = plot()
    for i = 1:npts
        valid = .~ismissing.(C[i,varname])
        if ~isinf(dtmax)
            valid = valid .& (C[i,"date_dt"] .<= dtmax)
        end
        p = plot!(C[i,"mid_date"][valid], C[i,varname][valid], seriestype = :scatter)
    end
    return p
end
