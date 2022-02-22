"""
lsqfit_interp(t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, t_i)

creates a continuous time series of velocity from the outputs of ITS_LIVE.lsqfit.


# Example
```julia
julia> v_i, v_i_err = lsqfit_interp(t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, t_i)
```

# Arguments
   - `t_fit:::Vector{DateTime}`: date of v_fit
   - `v_fit::Vector{Any}`: mean annual v
   - `amp_fit::Vector{Any}`: amplitude of seasonal cycle
   - `phase_fit::Vector{Any}`: phase of seasonal cycle
   - `v_fit_err::Vector{Any}`: error in mean annual v
   - `amp_fit_err::Vector{Any}`: error in amplitude of seasonal cycle
   - `t_i:::Vector{DateTime}`: dates to interpolate to

# Author
Chad A. Greene [original Matlab code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 1, 2022

Alex Gardner [Julia code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 17, 2022

"""
function lsqfit_interp(t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, t_i)
SplineOrder = 3;

# convert date to decimal years
t_fit = ITS_LIVE.decimalyear(t_fit)
t_i = ITS_LIVE.decimalyear(t_i)

valid = (.~ismissing.(t_fit)) .& (.~ismissing.(v_fit)) 
t_fit = convert.(Float64, t_fit[valid])
v_fit = convert.(Float64, v_fit[valid])
amp_fit = convert.(Float64, amp_fit[valid]) # not sure if the convert is needed
phase_fit = convert.(Float64, phase_fit[valid])

# interpolate discrete data using cubic splines (B-spline order k = 4)
itp = BSplineKit.interpolate(t_fit, v_fit, BSplineOrder(SplineOrder))
v_fit_i = itp.(t_i)  

if isa(amp_fit, Number)
   v_seas = sineval(hcat(amp_fit, ph),t_i); 
else
   # Convert from "polar" coordinates (where day-of-year phase is converted to radians) to cartesian coordinates so we can interpolate phase:  
   tmpx = amp_fit.*cos.(phase_fit*2*pi/365.25); # see "Algorithms" for https://www.mathworks.com/help/matlab/ref/pol2cart.html
   tmpy = amp_fit.*sin.(phase_fit*2*pi/365.25);
   
   # apply smooth interpoaltion in cartesian space
   itpx = BSplineKit.interpolate(t_fit, tmpx, BSplineOrder(SplineOrder))
   tmpxi = itpx.(t_i)
   itpy = BSplineKit.interpolate(t_fit,tmpy, BSplineOrder(SplineOrder))
   tmpyi = itpy.(t_i)  

   # Convert back to "polar" coordinates 
   tmp_phase_fit_i, amp_fit_i = atan.(tmpxi, tmpyi), hypot.(tmpxi,tmpyi)  # equivelent to Matlab's cart2pol(x,y)
   
   # phase_fit_i = tmp_phase_fit_i*365.25/(2*pi); 

   # create sine wave
   v_seas = amp_fit_i.*sin.(2*pi*t_i .+ tmp_phase_fit_i)
end

v_i = v_fit_i .+ v_seas; 

# compute error
itp = BSplineKit.interpolate(t_fit, vec(hypot.(v_fit_err, amp_fit_err)), BSplineOrder(SplineOrder))
v_i_err = itp.(t_i)  

return v_i, v_i_err
end