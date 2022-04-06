"""
    save2h5(C, outfilenamename)

[experimental] save named ITS_LIVE array (`C`) to an `outfilenamename` as an h5 file


# Example
```julia
julia> save2h5(C, outfilenamename)
```

# Author
Alex S. Gardner, JPL, Caltech.
"""


function save2h5(C, outfilename)
    # save data to a local HDF5
    # outfilename = "/Users/XXX/Desktop/ITS_LIVE_output.h5"
    h5open(outfilename, "w")
    varnames = names(C)[2];
    if ndims(C) == 2
        for pt = 1:size(C,1)
            for vname in varnames
                v = C[pt,vname]
                if v isa Vector{Union{Missing, Float64}}
                    v = Float64.(replace(v, missing => NaN))
                elseif vname == "satellite_img1" || vname == "satellite_img2"
                    id, sensorgroups = ItsLive.sensorgroup(v)
                    v = id;
                    sensor_info = "";
                    for i = 1:length(sensorgroups)
                        s = sensorgroups[i]
                        if i == length(sensorgroups)
                            sensor_info = sensor_info * repr(s["id"]) * "=" * s["name"]
                        else
                            sensor_info = sensor_info * repr(s["id"]) * "=" * s["name"] * ", " 
                        end
                    end
                    
                    write(outfilename, vname*"_info", sensor_info)
                end
                write(outfilename, vname * string(pt), v)
            end
        end
    else
        for vname in varnames
            v = C[vname]
            if v isa Vector{Union{Missing, Float64}}
                v = Float64.(replace(v, missing => NaN))
            elseif vname == "satellite_img1" || vname == "satellite_img2"
                id, sensorgroups = ItsLive.sensorgroup(v)
                v = id;
                sensor_info = "";
                for i = 1:length(sensorgroups)
                    s = sensorgroups[i]
                    if i == length(sensorgroups)
                        sensor_info = sensor_info * repr(s["id"]) * "=" * s["name"]
                    else
                        sensor_info = sensor_info * repr(s["id"]) * "=" * s["name"] * ", " 
                    end
                end
                
                write(outfilename, vname*"_info", sensor_info)
            end
            write(outfilename, vname * pt, v)
        end
    end
end