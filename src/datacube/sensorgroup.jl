"""
    id, sensorgroups= sensorgroup(sensor)

return the `sensor` group `id` and the corresponding `sensorgroups`


using Statistics

# Example
```julia
julia> id, sensorgroups = sensorgroup(sensor)
```

# Arguments
   - `sensor::Vector{Any}`: sensor list

# Author
Alex S. Gardner, JPL, Caltech.
"""
function sensorgroup(sensor)
    # specify groups of sensors to be filtered together
    sensorgroups = []


    # ensure that Sentinel 2 is the first index [has implications for dependencies]
    if ~(sensorgroups[1]["name"] == "Sentinel 2")
        error("Sentinel 2 is not the fist index of sensorgroups, check the sensorgroup function")
    end

    id = zeros(Int16, length(sensor))

    for sg = 1:length(sensorgroups)
        for k = 1:length(sensorgroups[sg]["sensors"])
            id[cmp.(sensorgroups[sg]["sensors"][k], sensor) .== 0] .= sg
        end
    end


    return id, sensorgroups
end