"""
    lsqfit(v,v_err,mid_date,date_dt,mad_thresh[optional],filt_iterations[optional])

error wighted model fit to discrete interval data 

using Statistics

# Example
```julia
julia> t_fit, v_fit, amp_fit, phase_fit, amp_fit_err, v_fit_err, fit_count, fit_outlier_frac = lsqfit(v,v_err,mid_date,date_dt,mad_thresh)
```

# Arguments
   - `v::Vector{Any}`: image-pair (discrete interval) velocity
   - `v_err::Vector{Any}`: image-pair (discrete interval) velocity error
   - `mid_date::Vector{DateTime}`: center date of image-pair [] 
   - `date_dt::Vector{Any}`: time seperation between image pairs [days]
   - `mad_thresh::Number`: optional input for MAD treshold for outlier rejection
   - `filt_iterations::Number`: optional input for number of iterations for MAD outlier rejection filter

# Author
Chad A. Greene [original Matlab code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 1, 2022

Alex Gardner [Julia code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 17, 2022

"""


function lsqfit(v, v_err, mid_date, date_dt, mad_thresh::Number = 5, filt_iterations::Number = 3)

#=
# add systimatic error based on level of co-registration
vy_error[stable_shift.==0] .= vy_error[stable_shift.==0] .+ 100
vx_error[stable_shift.==2] .= vx_error[stable_shift.==2] .+ 20
vy_error[stable_shift.==2] .= vy_error[stable_shift.==2] .+ 20
vx_error[stable_shift.==1] .= vx_error[stable_shift.==1] .+ 5
vy_error[stable_shift.==1] .= vy_error[stable_shift.==1] .+ 5
=#

outlier = ismissing.(v) 

t1 = mid_date .- (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))
t2 = mid_date .+ (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))

# Convert datenums to decimal years:
yr1 = ITS_LIVE.decimalyear(t1)
yr2 = ITS_LIVE.decimalyear(t2)

# dt in years:
dyr = yr2 .- yr1

# weights for velocities:
w_v = 1 ./ (v_err.^2)

# Weights (correspond to displacement error, not velocity error):
w_d = transpose(1. /( v_err .* dyr)) # Not squared because the p= line below would then have to include sqrt(w) on both accounts

d_obs = v.*dyr; # observed displacement in meters

#=
## Initial outlier detection and removal:
# I've tried this entire function with or without this section, and it
# doesn't seem to make much difference. It's actually fine to do
# the whole analysis without detrending at all, but I'm leaving this here
# for now anyways.
# #
# # # Define outliers as anything that's more than 10 standard deviations away from detrended vals.
# # outliers = abs(v) > 10*std(v);
# #
# # # Remove them!
# # yr = yr(~outliers,:);
# # dyr = dyr(~outliers);
# # d_obs = d_obs(~outliers);
# # w_d = w_d(~outliers);
# # w_v = w_v(~outliers);
=#


# apply an intitial w point running filter
valid = .!outlier
p = sortperm(mid_date[valid]);
w = 15;
vmed = FastRunningMedian.running_median((convert.(Float64, v[valid][p])),w)
resid = abs.(v[valid][p] - vmed);
sigma = Statistics.median(resid)*1.4826;

foo = @view outlier[valid];
foo[p[resid .> (mad_thresh*2*sigma)]] .= true # multiply threshold by 2 as this is a crude filter


## Make matrix of percentages of years corresponding to each displacement measurement
D, tD, M = ITS_LIVE.design_matrix(t1, t2, "interannual")
yr = ITS_LIVE.decimalyear(tD)

# Iterative mad filter []
for i = 1:filt_iterations
    valid = .!outlier
    
    # Solve for coefficients of each column in the Vandermonde:
    p = (w_d[valid].*D[valid,:]) \ (w_d[valid].*d_obs[valid]);

    ## Find and remove outliers    
    d_model = sum(broadcast(*,D[valid,:],transpose(p)),dims=2); # modeled displacements (m)
    
    d_resid = abs.(d_obs[valid] - d_model)./dyr[valid]; # devide by dt to avoid penalizing long dt [asg]
    
    d_sigma = Statistics.median(d_resid)*1.4826; # robust standard deviation of errors, using median absolute deviation

    # valid = vec(d_resid .<= (mad_thresh*d_sigma));
    outlier[valid] = vec(d_resid .> (mad_thresh*d_sigma))

    ## Remove no-data columns from M:
    #hasdata = vec(sum(M, dims = 1).>1);
    #yr = yr[hasdata];
    #M = M[:,hasdata];
end


D, tD, M = ITS_LIVE.design_matrix(t1, t2, "sinusoidal_interannual")
yr = ITS_LIVE.decimalyear(tD)

valid = .!outlier
fit_outlier_frac = sum(valid)./length(valid);

# Solve for coefficients of each column in the Vandermonde:
p = (w_d[valid].*D[valid,:]) \ (w_d[valid].*d_obs[valid]);

## Postprocess

# Convert coefficients to amplitude and phase of a single sinusoid:
Nyrs = length(yr);
amp_fit = hypot.(p[1:Nyrs],p[Nyrs+1:2*Nyrs]); # amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).

## THIS COULD BE SOURCE OF ERRORS AS MATLAB USED atan2
phase_rad = atan.(p[Nyrs+1:2*Nyrs],p[1:Nyrs]); # phase in radians
phase_fit = 365.25*(mod.(0.25 .- phase_rad/(2*pi),1)); # phase converted such that it reflects the day when value is maximized

# Goodness of fit:
d_model = sum(broadcast(*,D,transpose(p)),dims=2);

# A_err is the *velocity* (not displacement) error, which is the displacement error divided by the weighted mean dt:
amp_fit_err = Vector{Union{Float64,Missing}}(missing, size(amp_fit))

function stdw(x,w)
    μ = mean(x)
    s = sqrt(sum(w.*(x.-μ).^2)./sum(w))
    return s
end

for k = 1:Nyrs
    ind = (M[:,k] .> 0) .& valid;
    amp_fit_err[k] = stdw(d_obs[ind]-d_model[ind],w_d[ind]) ./ (sum(w_d[ind].*dyr[ind])./sum(w_d[ind])); # asg replaced call to wmean [!!! FOUND AND FIXED ERROR !!!!!!]
end

t_fit  = Dates.DateTime.(round.(Int,yr),7,1)
v_fit = p[2*Nyrs+1:end];

# Number of equivalent image pairs per year: (1 image pair equivalent means a full year of data. It takes about 23 16-day image pairs to make 1 year equivalent image pair.)
fit_count = sum(M[valid,:].>0, dims=1);

v_fit_err =  transpose(1 ./ sqrt.(sum(w_v[valid].*M[valid,:], dims=1)));

return t_fit, v_fit, amp_fit, phase_fit, amp_fit_err, v_fit_err, fit_count, fit_outlier_frac, outlier

end