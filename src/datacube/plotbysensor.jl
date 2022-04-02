"""
    p = plotbysensor(C, varname; dtmax::Number = Inf))

    plot its_live data (`C`) variable (`varname`) by sensor type 

# Example no inputs
```julia
julia> p = plotbysensor(x,y,sensor)
```

# Arguments
   - `C::Named Matrix{Any}`: ITS_LIVE named matrix [size(C,2) must = 1]
   - `varname::String`: variable name to plot
# Keyword Arguments
    - `dtmax::Number`: maximum time seperation between image pairs [days] to plot

# Author
Alex S. Gardner, JPL, Caltech.
"""
function plotbysensor(C, varname::String; dtmax::Number = Inf)

    if size(C,2) > 1
        error("plotbysensor() currently only accepts size(C,2) must = 1, try passing C[1,:] instead")
    end

    # determine sensor group ids
    valid = .~ismissing.(C[varname])
    
    if ~isinf(dtmax)
        valid = valid .& (C["date_dt"] .<= dtmax)
    end
    
    # find sensor groupings
    id, sensorgroups = ItsLive.sensorgroup(C["satellite_img1"][valid])

    # unique sensors
    uid = unique(id)

    p = plot()
    for i = 1:length(uid)
        ind = id .== uid[i]
        plot!(C["mid_date"][valid][ind], C[varname][valid][ind], seriestype = :scatter,  label = sensorgroups[uid[i]]["name"])
    end
    display(p)
    return p
end