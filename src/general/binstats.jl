"""
    binstats(x, y; binedges::Vector{Float64} = [0.0], dx::Number = 0, method::String = "mean")

    return statistics of `x` central value and `x` spread according to `method` ["mean" = default]
    argument on values binned by `y`.

# Example no inputs
```julia
julia> yb, ybstd, countb, bincenters = binstats(x, y; dx = 1)
```

# Arguments
   - `x::Vector{Float64}`: x data [typically DateTime]
   - `y::Vector{Float64}`: data to be binned
   - `dt::Number`: width of x bins [dt *or* binedges must be provided]
   - `binedges::Vector{Float64}`: x bin edges [dt *or* binedges must be provided]
   - `method::String`: method for deter
   - `skipspread::Bool`: flag for skipping spread calculations

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 25, 2022
"""

function binstats(x, y; binedges::Vector{Float64} = [0.0], dx::Number = 0, method::String = "mean", skipspread::Bool = false)

    if (dx != 0) & (length(binedges) == 1)
        # make bin edges
        binedges = ((floor(minimum(x)/(dx/2))-1)*(dx/2)):dx:((ceil(maximum(x)/(dx/2)) + 1)*(dx/2))
    elseif (length(binedges) > 1) & (dx == 0)
        # bin edges are provided
    elseif (dx == 0) & (length(binedges) == 1)
        error("dx *or* binedges must be provided")
    else
        error("only dx *or* binedges should be provided, not both")
    end

    yb = fill(NaN, length(binedges)-1)
    if !skipspread
        ybstd = copy(yb)
    end
    countb = zeros(Int, length(yb))

    valid = .~ismissing.(y)
    x = x[valid]
    y = y[valid]

    # sort data
    sind = sortperm(x)
    x = x[sind]
    y = y[sind]

    # find bin index
    ind0 = searchsortedfirst.(Ref(binedges[1:end]), x)
    
    for i = 1:(length(binedges)-1) 
        ind = searchsorted(ind0, i+1) # the +1 is so that data <= binedges[1] is not included
        if !isempty(ind)
            countb[i] = sum(ind.stop-ind.start+1)
            if method == "mean"

                yb[i] = mean(y[ind])
                # calculate standard deviation
                if !skipspread
                    ybstd[i] = sqrt(sum((y[ind] .- yb[i]).^2)/countb[i])
                end
            elseif method == "median"
                yb[i] = Statistics.median(y[ind])
                # MAD scaled to match standard deviation
                if !skipspread
                    ybstd[i] = Statistics.median(y[ind] .- yb[i])*1.4826
                end
            else
                error("method not recognized")     
            end  
        end
    end

    # bin centers
    bincenters = (binedges[1:end-1] + binedges[2:end])./2
    
    if !skipspread
        return yb, ybstd, countb, bincenters
    else
        return yb, countb, bincenters
    end
end