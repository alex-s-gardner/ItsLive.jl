"""
    running_mean(v, w)

return the running mean of `v` using kernel width `w`
using StatiStatistics

# Example
```julia
julia> running_mean(v, w)
```

# Arguments
   - `v`: input data
   - `w::Number`: width of running mean

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 18, 2022
"""

function running_mean(v, w::Number) 
    out = [i < w ? Statistics.mean(v[begin:i]) : Statistics.mean(v[i-w+1:i]) for i in 1:length(v)]
    return out
end