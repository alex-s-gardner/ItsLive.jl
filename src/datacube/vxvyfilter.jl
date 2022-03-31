"""
    outlier, dtmax, sensorgroups = vxvyfilter(vx,vy,dt; sensor)

identify 'outlier's for which the `vx` or `vy` distibution changes for longer `dt`. 'dtmax' identified for each 'sensorgroups'

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
Alex S. Gardner, JPL, Caltech.
"""
function vxvyfilter(vx,vy,dt; sensor::Vector = ["none"])
    # filter bins
    binedges = [0, 16, 32, 64, 128, 256, 1E10]

    # initialize output
    outlier = falses(size(vx))

    # project vx and vy onto the median flow vector for dt <= 32
    ind = (dt .<= 16) .& .~ismissing.(vx)
    vx0 = Statistics.median(vx[ind])
    vy0 = Statistics.median(vy[ind])

    # check that vx0 and vy0 are not both equal to zero... if so set vx0 = 1;
    # this is equivelent to picking an arbitrary direction on which to project 
    #component velocities and avoids devide by zero  
    if (vx0 == 0) & (vy0 == 0)
        vx0 = 1;
        vx0 = 1;
    end

    v0 = sqrt(vx0.^2 + vy0.^2);
    uv = vcat(vx0/v0, vy0/v0) # unit flow vector
    vp = hcat(vx, vy) * uv # flow acceleration in direction of unit flow vector

    # if sensor variable is included then seperate filtering by sensor
    if sensor[1] == "none"
        # initialize output
        dtmax = Vector{Union{Missing, Float64}}(missing, 1)
        
        ind = .~ismissing.(vp)

        # find the maximum dt that is not significantly different from the minimum dt bin
        vpdtmax = dtfilter(vp[ind], dt[ind], binedges)

        # find the minimum acceptable dt threshold
        dtmax[1] = vpdtmax

        if dtmax[1] > 20E3
            # no data needs to be masked
        else
            # replace data that exceed dt threshold with missings 
            outlier = ((dt .> dtmax[1]) .& ind) .| outlier
        end

        sensorgroups = Vector{Vector{String}}(undef,1)
        sensorgroups[1] = ["none"]

    else
        # specify groups of sensors to be filtered together
        id, sensorgroups = sensorgroup(sensor)

        # initialize output
        dtmax = Vector{Union{Missing, Float64}}(missing, size(sensorgroups))

        for sg = 1:length(sensorgroups)
            sgind = id .== sg
            
            if any(sgind)
                #find the maximum dt that is not significantly different from the minimum dt bin
                vpdtmax = dtfilter(vp[sgind], dt[sgind], binedges)

                # find the minimum acceptable dt threshold
                dtmax[sg] = vpdtmax

                if dtmax[sg] > 20E3
                    # no data needs to be masked
                else
                    # replace data that exceed dt threshold with missings 
                    valind = (.~ismissing.(vp)) .& sgind
                    outlier = ((dt .> dtmax[sg]) .& valind) .| outlier 
                end
            end
        end
    end
    return outlier, dtmax, sensorgroups
    
end
