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
    push!(sensorgroups, Dict([("id", 1), ("name", "Sentinel 2"), ("sensors", ["2A", "2B"])]))
    push!(sensorgroups, Dict([("id", 2), ("name", "Landsat 8/9"), ("sensors",["8.0", "9.0"])]))
    push!(sensorgroups, Dict([("id", 3), ("name", "Sentinel 1"), ("sensors",["1A", "1B"])]))
    push!(sensorgroups, Dict([("id", 4), ("name", "Landsat 7"), ("sensors",["7.0"])]))
    push!(sensorgroups, Dict([("id", 5), ("name", "Landsat 4/5"), ("sensors",["4.0", "5.0"])]))

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