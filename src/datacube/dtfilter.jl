"""
    dtfilter(x,dt,binedges)

returns the maximum dt for which the distribution of x shows no statistical difference from the distribution of x in the minimum dt bin

This filter is needed to identify longer dts that exhibit "skipping" or "locking" behavior in feature tracking estimates of surface flow. This happens when the surface texture provides a lesser match than to stationary features, due to long time separation between repeat images, such as ice falls and curved medial moraines.

using Statistics

# Example
```julia
julia> dtfilter(vx,dt,binedge)
```

# Arguments
   - `x::Vector{Any}`: value to be filtered as a function of dt (typically velocity)
   - `dt::Vector{Any}`: time seperation between image pairs [days]
   - `binedges::Vector{Float64}`: *optional* edges of dt bins into which vx and vy will be grouped and compared
   - `dtbin_mad_thresh::Number`: used to determine in dt means are significantly different

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 10, 2022
"""

function dtfilter(x, dt , binedges::Vector{Float64} = [0, 32, 64, 128, 256, 1E10], dtbin_mad_thresh::Number = 2)


    ## define internal functions 
    # define a median absolut difference function
    function medmad(x) 
        if all(ismissing, x)
            medx = 0
            madx = 0
        else
            medx = round(Int,Statistics.median(skipmissing(x)));
            madx = round(Int,Statistics.median(skipmissing(abs.(x .- medx))));
        end
        return [medx, madx]
    end

    # find last valid bin or return zero
    function findlastorzero(x::BitMatrix)
        if any(x)
            x = findlast(x)
        else
            x = 0
        end
        return x
    end

    # round dt to an even integer
    dt = round.(Int, dt)

    # sort dt and x bu decending dt
    dt_sortperm = sortperm(dt);
    dt = dt[dt_sortperm]
    x = x[dt_sortperm]

    # find the last sorted value that is < 
    bin_ind = searchsortedlast.(Ref(dt), binedges) 
    bin_ind[1] = 1;

    binMad = zeros(Int16, length(binedges)-1)
    binMed = zeros(Int16, length(binedges)-1)

    size(binMad)
    for i in range(1, length(bin_ind)-1)
            foo = medmad(x[bin_ind[i]:bin_ind[i+1]]);
            binMed[i] = foo[1]
            binMad[i] = foo[2]
    end

    # check if populations overlap (use first, smallest dt, bin as reference)
    minBound = binMed - (binMad * dtbin_mad_thresh * 1.4826);
    maxBound = binMed + (binMad * dtbin_mad_thresh * 1.4826);

    # check if distributions overlap
    exclude = (minBound .> maxBound[1]) .| (maxBound .< minBound[1]);

    if !any(exclude)
        dtmax = (2^15 - 1)
    else
        dtmax = findfirst(exclude)
        dtmax = binedges[dtmax]
    end
    return dtmax
end