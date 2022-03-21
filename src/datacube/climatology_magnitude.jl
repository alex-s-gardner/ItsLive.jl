"""
    climatology_magnitude(vx0, vy0, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err, vx_phase, vy_phase)

returns the mean, trend, seasonal amplitude, error in seasonal amplitude, and seasonal phase from component values projected on the unit flow  vector defined by vx0 and vy0

# Example no inputs
```julia
julia> v, dv_dt, v_amp, v_amp_err, v_phase = climatology_magnitude(vx0, vy0, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err, vx_phase, vy_phase)
```

# Arguments
   - `vx0::Number`: mean flow in x direction
   - `vy0::Number`: mean flow in y direction
   - `dvx_dt::Number`: trend in flow in x direction
   - `dvy_dt::Number`: trend in flow in y direction
   - `vx_amp::Number`: seasonal amplitude in x direction
   - `vy_amp::Number`: seasonal amplitude in y direction 
   - `vx_amp_err::Number`: error in seasonal amplitude in x direction
   - `vy_amp_err::Number`: error in seasonal amplitude in y direction 
   - `vx_phase::Number`: seasonal phase in x direction [day of maximum flow]
   - `vy_phase::Number`: seasonal phase in y direction [day of maximum flow]


# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 20, 2022
"""

function climatology_magnitude(vx0, vy0, dvx_dt, dvy_dt, vx_amp, vy_amp, vx_amp_err, vy_amp_err, vx_phase, vy_phase)
    
    # solve for velcity magnitude and acceleration [do this using vx and vy as to not bias the result due to the Rician distribution of v]
    v = sqrt(vx0^2 + vy0^2) # velocity magnitude
    uv = vcat(vx0/v, vy0/v) # unit flow vector
    dv_dt = hcat(dvx_dt, dvy_dt) * uv # flow acceleration in direction of unit flow vector
    dv_dt = dv_dt[1] # conver from vector to number
    v_amp_err = abs.(hcat(vx_amp_err, vy_amp_err)) * abs.(uv) # flow acceleration in direction of unit flow vector, take absolute values
    v_amp_err = v_amp_err[1] # conver from vector to number 
    # solve for amplitude and phase in unit flow dirction
    if isa(vx_amp, Number)
        t0 = 0:0.1:1
        vx_sin = vx_amp .* sin.((t0 .+ (-vx_phase/365.25 + 0.25))*2*pi);  # must convert phase to fraction of a year and adjust from peak to phase
        vy_sin = vy_amp .* sin.((t0 .+ (-vy_phase/365.25 + 0.25))*2*pi);  # must convert phase to fraction of a year and adjust from peak to phase
        v_sin = hcat(vx_sin, vy_sin) * uv # seasonality in direction of unit flow vector
        a1, a2 = hcat(cos.(t0 .* (2*pi)) ,  sin.(t0 .* (2*pi)))\ v_sin
        v_amp = hypot(a1, a2); # amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).
        phase_rad = atan.(a1,a2); # phase in radians
        v_phase = 365.25*(mod.(0.25 .- phase_rad/(2*pi),1)); # phase converted such that it reflects the day when value is maximized
    else
        error("amplitude and phase in unit flow dirction not yet implimented for multiple phase and amplitude")
    end

    return v, dv_dt, v_amp, v_amp_err, v_phase
end