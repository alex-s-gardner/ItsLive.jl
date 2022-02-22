"""
    lsqfit_interp(t_int,v_int,amp,phase,ti,v_int_err,amp_err)

creates a continuous time series of velocity from the outputs of lsqfit.


# Example
```julia
julia> vi, v_int_i = lsqfit_interp(t_int,v_int,amp,phase,ti,v_int_err,amp_err)
```

# Arguments
   - `t_int:::Vector{DateTime}`: date of v_int
   - `v_int::Vector{Any}`: mean annual v
   - `amp::Vector{Any}`: amplitude of seasonal cycle
   - `phase::Vector{Any}`: phase of seasonal cycle
   - `ti:::Vector{DateTime}`: dates to interpolate to
   - `v_intIerr::Vector{Any}`: error in mean annual v
   - `amp_err::Vector{Any}`: error in amplitude of seasonal cycle

# Author
Chad A. Greene [original Matlab code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 1, 2022

Alex Gardner [Julia code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 17, 2022

"""
function lsqfit_interp(t_int,v_int,amp,phase,ti,v_int_err,amp_err)
SplineOrder = 3;

# convert date to decimal years
t_int = ITS_LIVE.decimalyear(t_int)
ti = ITS_LIVE.decimalyear(ti)

valid = (.~ismissing.(t_int)) .& (.~ismissing.(v_int)) 
t_int = convert.(Float64, t_int[valid])
v_int = convert.(Float64, v_int[valid])
amp = convert.(Float64, amp[valid]) # not sure if the convert is needed
phase = convert.(Float64, phase[valid])

# interpolate discrete data using cubic splines (B-spline order k = 4)
itp = BSplineKit.interpolate(t_int, v_int, BSplineOrder(SplineOrder))
v_int_i = itp.(ti)  

if isa(amp, Number)
   v_seas = sineval(hcat(amp, ph),ti); 
else
   # Convert from "polar" coordinates (where day-of-year phase is converted to radians) to cartesian coordinates so we can interpolate phase:  
   tmpx = amp.*cos.(phase*2*pi/365.25); # see "Algorithms" for https://www.mathworks.com/help/matlab/ref/pol2cart.html
   tmpy = amp.*sin.(phase*2*pi/365.25);
   
   # apply smooth interpoaltion in cartesian space
   itpx = BSplineKit.interpolate(t_int, tmpx, BSplineOrder(SplineOrder))
   tmpxi = itpx.(ti)
   itpy = BSplineKit.interpolate(t_int,tmpy, BSplineOrder(SplineOrder))
   tmpyi = itpy.(ti)  

   # Convert back to "polar" coordinates 
   tmp_phase_i, amp_i = atan.(tmpxi, tmpyi), hypot.(tmpxi,tmpyi)  # equivelent to Matlab's cart2pol(x,y)
   
   # phase_i = tmp_phase_i*365.25/(2*pi); 

   # create sine wave
   v_seas = amp_i.*sin.(2*pi*ti .+ tmp_phase_i)
end

vi = v_int_i .+ v_seas; 

# compute error
itp = BSplineKit.interpolate(t_int, vec(hypot.(v_int_err, amp_err)), BSplineOrder(SplineOrder))
v_int_i = itp.(ti)  

return vi, v_int_i
end