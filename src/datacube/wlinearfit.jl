"""
    offset, slope, error = wliearfit(t, v, v_err, datetime0)

return the `offset`, `slope`, and `error` for a weighted linear fit to `v` with an intercept of `datetime0`

# Example no inputs
```julia
julia> offset, slope, error = wliearfit(t, v, v_err, datetime0)
```

# Arguments
   - `t::Vector{DateTime}`: date of input estimates
   - `v::Vector{Float64}`: estimates
   - `v_err::Vector{Float64}`: estimate errors
   - `datetime0::DateTime)`: model intercept

# Author
Alex S. Gardner, JPL, Caltech.
"""
function wlinearfit(t::Vector{DateTime}, v::Vector{Float64} , v_err::Vector{Float64} , datetime0::DateTime)

    yr = ItsLive.decimalyear(t)
    yr0 = ItsLive.decimalyear(datetime0)

    yr = yr .- yr0;

    # weights for velocities:
    w_v = 1 ./ (v_err.^2)

    # create design matrix
    D = ones((length(yr),2))
    D[:,2] = yr;

    # Solve for coefficients of each column in the design matrix
    valid = .!ismissing.(v)
    offset, slope = (w_v[valid].*D[valid,:]) \ (w_v[valid].*v[valid]);

    # fit error
    error = sqrt(sum(v_err[valid].^2))/(sum(valid)-1)

    return offset, slope, error
end 
