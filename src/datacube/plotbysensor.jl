"""
    plotbysensor(x,y,sensor)

    plot ITS_LIVE data `x` and `y` colored by `sensor`

# Example no inputs
```julia
julia> p = plotbysensor(x,y,sensor)
```

# Arguments
   - `x::Vector{Float64}`: x data [typically DateTime]
   - `y::Vector{Float64}`: y data
   - `sensor::Vector{Any}`: sensor

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 27, 2022
"""

function plotbysensor(x,y,sensor)

    # determine sensor group ids
    valid = .~ismissing.(x)
    if any(.~valid)
        sensor = sensor[valid]
        x = x[valid]
        y = y[valid]
    end

    # find sensor groupings
    id, sensorgroups = ItsLive.sensorgroup(sensor)

    # unique sensors
    uid = unique(id)

    p = plotly()
    plot()
    for i = 1:length(uid)
        ind = id .== uid[i]
        p = plot!(x[ind], y[ind], seriestype = :scatter,  label = sensorgroups[uid[i]]["name"])
    
    end
    display(p)
    return p
end