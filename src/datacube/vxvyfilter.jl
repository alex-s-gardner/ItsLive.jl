"""
    vxvyfilter(vx,vy,dt,sensor)

remove data for which the vx or vy distibution changes for longer dts

This filter is needed to identify longer dts that exhibit "skipping" or "locking" behavior in feature tracking estimates of surface flow. This happens when the surface texture provides a lesser match than to stationary features, due to long time separation between repeat images, such as ice falls and curved medial moraines.

using Statistics

# Example
```julia
julia> outlier, dtmax, sensorgroups = vxvyfilter(vx,vy,dt,sensor)
```

# Arguments
   - `vx::Vector{Any}`: x component velocity
   - `vy::Vector{Any}`: y component velocity
   - `dt::Vector{Any}`: time seperation between image pairs [days]
   - `sensor::Vector{Any}`: list of image sensors for image-pairs... used to group results

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 10, 2022
"""

function vxvyfilter(vx,vy,dt,sensor::Vector{Any})

    # specify groups of sensors to be filtered together
    sensorgroups = Vector{Vector{String}}(undef,1)
    sensorgroups[1] = ["8.", "9."]
    push!(sensorgroups, ["2A", "2B"])
    push!(sensorgroups, ["1A", "1B"])
    push!(sensorgroups, ["7."])
    push!(sensorgroups, ["4.", "5."])

    # filter bins
    binedges = [0, 32, 64, 128, 256, 1E10]

    # initialize output
    dtmax = Vector{Union{Missing, Float64}}(missing, size(sensorgroups))
    outlier = falses(size(vx))
    for sg = 1:length(sensorgroups)

        sgind = falses(size(sensor))
        for k = 1:length(sensorgroups[sg])
            sgind = sgind .| (sensor .== sensorgroups[sg][k])
        end
        
        if any(sgind)

            #find the maximum dt that is not significantly different from the minimum dt bin
            vxdtmax = dtfilter(vx[sgind], dt[sgind], binedges)
            vydtmax = dtfilter(vy[sgind], dt[sgind], binedges)

            # find the minimum acceptable dt threshold
            dtmax[sg] = min(vxdtmax, vydtmax)

            if dtmax[sg] > 20E3
                # no data needs to be masked
            else
                # replace data that exceed dt threshold with missings 
                valind = (.~ismissing.(vx)) .& sgind
                outlier = ((dt .> dtmax[sg]) .& valind) .| outlier 
            end
        end
    end
    return outlier, dtmax, sensorgroups
end
