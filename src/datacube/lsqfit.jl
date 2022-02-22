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


function lsqfit(v, v_err, mid_date, date_dt, mad_thresh::Number = 3, filt_iterations::Number = 3;)

#=
# add systimatic error based on level of co-registration
vy_error[stable_shift.==0] .= vy_error[stable_shift.==0] .+ 100
vx_error[stable_shift.==2] .= vx_error[stable_shift.==2] .+ 20
vy_error[stable_shift.==2] .= vy_error[stable_shift.==2] .+ 20
vx_error[stable_shift.==1] .= vx_error[stable_shift.==1] .+ 5
vy_error[stable_shift.==1] .= vy_error[stable_shift.==1] .+ 5
=#


valid = .~ismissing.(v) 

t1 = mid_date .- (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))
t2 = mid_date .+ (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))
t1 = t1[valid];
t2 = t2[valid];
v = v[valid]
v_err = v_err[valid]


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

## Make matrix of percentages of years corresponding to each displacement measurement

y1 = floor(minimum(yr1)):floor(maximum(yr2));

M = zeros(length(dyr),length(y1));

# Loop through each year:
for k = 1:length(y1)
    
    # Set all measurements that begin before the first day of the year and end after the last
    # day of the year to 1:
    ind = (yr1.<=y1[k]) .& (yr2.>=(y1[k]+1));
    M[ind,k] .= 1;
    
    # Within year:
    ind = (yr1.>=y1[k]) .& (yr2.<(y1[k]+1));
    M[ind,k] = dyr[ind];
    
    # Started before the beginning of the year and ends during the year:
    ind = (yr1.<y1[k]) .& (yr2.>=y1[k]) .& (yr2.<(y1[k].+1));
    M[ind,k] = yr2[ind] .- y1[k];
    
    # Started during the year and ends the next year:
    ind = (yr1.>=y1[k]) .& (yr1.<(y1[k].+1)) .& (yr2.>=(y1[k]+1));
    M[ind,k] = (y1[k]+1) .- yr1[ind];
end

hasdata = vec(sum(M,dims=1).>0);
y1 = y1[hasdata];
M = M[:,hasdata];

# Iterative mad filter []
totalnum = length(yr1);
for i = 1:filt_iterations

    # Displacement Vandermonde matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
    N = ((sin.(2*pi*yr2) - sin.(2*pi*yr1))./(2*pi));
    N = mapslices(row->(N).*row,float(M.>0), dims=1);

    P = (cos.(2*pi*yr1) - cos.(2*pi*yr2)) ./ (2*pi);
    P = mapslices(row->(P).*row,float(M.>0), dims=1);

    D = hcat(P, N, M);

    # Solve for coefficients of each column in the Vandermonde:
    p = (w_d.*D) \ (w_d.*d_obs);

    #= 
    # Displacement Vandermonde matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
    D = hcat((cos.(2*pi*yr1) - cos.(2*pi*yr2)) ./ (2*pi), (sin.(2*pi*yr2) - sin.(2*pi*yr1))./(2*pi), M, ones(size(dyr)));
    
    # Solve for coefficients of each column in the Vandermonde:
    p = (w_d.*D) \ (w_d.*d_obs);
    =#

    ## Find and remove outliers
    
    d_model = sum(broadcast(*,D,transpose(p)),dims=2); # modeled displacements (m)
    
    d_resid = abs.(d_obs - d_model)./dyr; # devide by dt to avoid penalizing long dt [asg]
    
    d_sigma = Statistics.median(d_resid)*1.4826; # robust standard deviation of errors, using median absolute deviation
    valid = vec(d_resid .<= (mad_thresh*d_sigma));

    # Remove outliers:
    yr1 = yr1[valid];
    yr2 = yr2[valid];
    dyr = dyr[valid];
    d_obs = d_obs[valid];
    w_d = w_d[valid];
    w_v = w_v[valid];
    M = M[valid,:];
    
    # Remove no-data columns from M:
    hasdata = vec(sum(M, dims = 1).>1);
    y1 = y1[hasdata];
    M = M[:,hasdata];
end

fit_outlier_frac = (totalnum-length(yr1))./totalnum;

## Second iteration

# Displacement Vandermonde matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
N = ((sin.(2*pi*yr2) - sin.(2*pi*yr1))./(2*pi));
N = mapslices(row->(N).*row,float(M.>0), dims=1);

P = (cos.(2*pi*yr1) - cos.(2*pi*yr2)) ./ (2*pi);
P = mapslices(row->(P).*row,float(M.>0), dims=1);

D = hcat(P, N, M);

# Solve for coefficients of each column in the Vandermonde:
p = (w_d.*D) \ (w_d.*d_obs);

## Postprocess

# Convert coefficients to amplitude and phase of a single sinusoid:
Nyrs = length(y1);
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
    ind = M[:,k] .> 0;
    amp_fit_err[k] = stdw(d_obs[ind]-d_model[ind],w_d[ind]) ./ (sum(w_d[ind].*dyr[ind])./sum(w_d[ind])); # asg replaced call to wmean [!!! FOUND AND FIXED ERROR !!!!!!]
end

t_fit  = Dates.DateTime.(round.(Int,y1),7,1)
v_fit = p[2*Nyrs+1:end];

# Number of equivalent image pairs per year: (1 image pair equivalent means a full year of data. It takes about 23 16-day image pairs to make 1 year equivalent image pair.)
fit_count = sum(M.>0);

v_fit_err =  transpose(1 ./ sqrt.(sum(w_v.*M, dims=1)));


return t_fit, v_fit, amp_fit, phase_fit, amp_fit_err, v_fit_err, fit_count, fit_outlier_frac

end