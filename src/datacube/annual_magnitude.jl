"""
    v_fit, v_fit_err, v_fit_count, v_fit_outlier_frac  = annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, vy_fit_err, vx_fit_count, vy_fit_count, vx_fit_outlier_frac, vy_fit_outlier_frac)

returns the annual mean, error, count, and outlier fraction from component values projected on the unit flow  vector defined by vx0 and vy0

# Example no inputs
```julia
julia>  annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, vy_fit_err, vx_fit_count, vx_fit_count, vx_fit_outlier_frac, vy_fit_outlier_frac)
```

# Arguments
    - `vx0::Number`: mean flow in x direction
    - `vy0::Number`: mean flow in y direction
    - `vx_fit::Vector`: annual mean flow in x direction
    - `vy_fit::Vector`: annual mean flow in y direction
    - `vx_fit_err::Vector`: error in annual mean flow in x direction
    - `vy_fit_err::Vector`: error in annual mean flow in y direction
    - `vx_fit_count::Vector`: number of values used to determine annual mean flow in x direction
    - `vy_fit_count::Vector`: number of values used to determine annual mean flow in y direction
    - `vx_fit_outlier_frac::Vector`: fraction of data identified as outliers and removed when calculating annual mean flow in x direction
    - `vy_fit_outlier_frac::Vector`: fraction of data identified as outliers and removed when calculating annual mean flow in y direction

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 20, 2022
"""

function annual_magnitude(vx0, vy0, vx_fit, vy_fit, vx_fit_err, vy_fit_err, vx_fit_count, vy_fit_count, vx_fit_outlier_frac, vy_fit_outlier_frac)
    
    # solve for velcity magnitude 
    v_fit = sqrt.(vx_fit.^2 + vy_fit.^2) # velocity magnitude
    uv = vcat(vx0/v_fit, vy0/v_fit) # unit flow vector
    v_fit_err = abs.(hcat(vx_fit_err, vy_fit_err)) * abs.(uv) # flow acceleration in direction of unit flow vector, take absolute values
    v_fit_count = ceil.((vx_fit_count .+ vy_fit_count) ./ 2)
    v_fit_outlier_frac = (vx_fit_outlier_frac .+ vy_fit_outlier_frac) ./ 2

    return v_fit, v_fit_err, v_fit_count, v_fit_outlier_frac
end