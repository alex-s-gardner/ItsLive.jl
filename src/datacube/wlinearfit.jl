"""
    wliearfit(t, v, v_err, datetime0)

returns the offset, slope, and se for a weighted linear fit to v with an intercept of datetime0

# Example no inputs
```julia
julia> offset, slope, se = wliearfit(t, v, v_err, datetime0)
```

# Arguments
   - `t::Vector{DateTime}`: date of input estimates
   - `v::Vector{Float64}`: estimates
   - `v_err::Vector{Float64}`: estimate errors
   - `datetime0::DateTime)`: model intercept

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 17, 2022
"""


function wlinearfit(t::Vector{DateTime}, v::Vector{Float64} , v_err::Vector{Float64} , datetime0::DateTime)

    yr = ITS_LIVE.decimalyear(t)
    yr0 = ITS_LIVE.decimalyear(datetime0)

    yr = yr .- yr0;

    # weights for velocities:
    w_v = 1 ./ (v_err.^2)

    # create design matrix
    D = ones((length(yr),2))
    D[:,2] = yr;

    # Solve for coefficients of each column in the design matrix
    valid = .!ismissing.(v)
    offset, slope = (w_v[valid].*D[valid,:]) \ (w_v[valid].*v[valid]);

    # RMSE from fit
    res = v .- (yr.*slope .+ offset);
    se = sqrt(sum(res.^2) ./ (sum(valid)-1))

    return offset, slope, se
end 
