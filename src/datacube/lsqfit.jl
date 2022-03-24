"""
    lsqfit(v,v_err,mid_date,date_dt; mad_thresh[optional], iterations[optional], model[optional])

error wighted model fit to discrete interval data 

using Statistics

# Example
```julia
julia> t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, fit_count, fit_outlier_frac, outlier = 
    lsqfit(v, v_err, mid_date, date_dt; mad_thresh, iterations, model)
```

# Arguments
   - `v::Vector{Any}`: image-pair (discrete interval) velocity
   - `v_err::Vector{Any}`: image-pair (discrete interval) velocity error
   - `mid_date::Vector{DateTime}`: center date of image-pair [] 
   - `date_dt::Vector{Any}`: time seperation between image pairs [days]
   - `mad_thresh::Number`: optional key word argument for MAD treshold for outlier rejection
   - `iterations::Number`: optional key word argument for number of iterations for lsqfit outlier rejection
   - `model::String`: optional key word argument giving name of lsq model. Valid options are "sinusoidal_interannual", "sinusoidal", "interannual".

# Author
Alex Gardner [Julia code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 17, 2022

Chad A. Greene [original Matlab inspriation code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 1, 2022
"""


function lsqfit(v, v_err, mid_date, date_dt; mad_thresh::Number = 6, iterations::Number = 1, model::String = "interannual")

outlier = ismissing.(v) 

t1 = mid_date .- (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))
t2 = mid_date .+ (Dates.Second.(round.(Int64,date_dt .* 86400 ./2)))

# Convert datenums to decimal years:
yr1 = ItsLive.decimalyear(t1)
yr2 = ItsLive.decimalyear(t2)
yr = floor(minimum(yr1)):floor(maximum(yr2));

# dt in years:
dyr = yr2 .- yr1

# weights for velocities:
w_v = 1 ./ (v_err.^2)

# Weights (correspond to displacement error, not velocity error):
w_d = transpose(1. /( v_err .* dyr)) # Not squared because the p= line below would then have to include sqrt(w) on both accounts
w_d = ones(size(w_d))

## pre filter data
valid = .!outlier
p = sortperm(mid_date[valid]);

# moving window MAD filter seems much more robust
w = 15;
resid = zeros(size(p))
vmed = FastRunningMedian.running_median((convert.(Float64, v[valid][p])),w)
resid[p] = abs.(v[valid][p] - vmed);

# this is the original polynomial fit that seems to have issues when there is strong seasonality
# yr_mid = (yr1+yr2)./2
# ply = Polynomials.fit(yr_mid[valid][p],collect(skipmissing(v[valid][p])), 2)
# resid = abs.(v[valid] .- ply.(yr_mid[valid]))

# filter data ouside of threshold
sigma = Statistics.median(resid)*1.4826;
foo = view(outlier, valid)
foo[resid .> (mad_thresh*2*sigma)] .= true # multiply threshold by 2 as this is a crude filter

# uncommnet following three lined to plot residules and outliers
# plot(mid_date[valid], resid, seriestype = :scatter, mc = :gray, labels="residuals")
# pp = plot!(mid_date[outlier[valid]], resid[outlier[valid]], seriestype = :scatter, mc = :red, labels="outliers")
# display(pp)

# Iterative mad filter []

## Make matrix of percentages of years corresponding to each displacement measurement
D, tD, M = ItsLive.design_matrix(t1, t2, model)
validYears = vec(sum(M[valid,:],dims=1) .> 1)
if model == "sinusoidal_interannual"
    validD = repeat(validYears,3)
elseif model == "sinusoidal"
    validD = vcat(BitVector([true, true]), validYears);
else
    validD = validYears
end

D = D[:,validD];
M = M[:,validYears];
yr = yr[validYears];

d_obs = v.*dyr; # observed displacement in meters
for i = 1:iterations

    # Solve for coefficients of each column in the Vandermonde:
    p = (w_d[valid].*D[valid,:]) \ (w_d[valid].*d_obs[valid]);

    if i < iterations
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
end

valid = .!outlier
fit_outlier_frac = 1 -(sum(valid)./length(valid));

# Goodness of fit:
d_model = sum(broadcast(*,D,transpose(p)),dims=2);

function stdw(x,w)
    μ = mean(x)
    s = sqrt(sum(w.*(x.-μ).^2)./sum(w))
    return s
end

# Convert coefficients to amplitude and phase of a single sinusoid:
Nyrs = length(yr);
if model == "sinusoidal_interannual"
    amp_fit = hypot.(p[1:Nyrs],p[Nyrs+1:2*Nyrs]); # amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).
    phase_rad = atan.(p[Nyrs+1:2*Nyrs],p[1:Nyrs]); # phase in radians
    phase_fit = 365.25*(mod.(0.25 .- phase_rad/(2*pi),1)); # phase converted such that it reflects the day when value is maximized

    # A_err is the *velocity* (not displacement) error, which is the displacement error divided by the weighted mean dt:
    amp_fit_err = Vector{Union{Float64,Missing}}(missing, size(amp_fit))
    for k = 1:Nyrs
        ind = (M[:,k] .> 0) .& valid;
        amp_fit_err[k] = stdw(d_obs[ind]-d_model[ind],w_d[ind]) ./ (sum(w_d[ind].*dyr[ind])./sum(w_d[ind])); # asg replaced call to wmean [!!! FOUND AND FIXED ERROR !!!!!!]
    end

    t_fit  = Dates.DateTime.(round.(Int,yr),7,1)
    v_fit = p[2*Nyrs+1:end];

elseif model == "sinusoidal"
    println(p)
    amp_fit = hypot.(p[1],p[2]); # amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).
    phase_rad = atan.(p[2],p[1]); # phase in radians
    phase_fit = 365.25*(mod.(0.25 .- phase_rad/(2*pi),1)); # phase converted such that it reflects the day when value is maximized

    # A_err is the *velocity* (not displacement) error, which is the displacement error divided by the weighted mean dt:
    amp_fit_err = Vector{Union{Float64,Missing}}(missing, 1)
    ind = valid;
    amp_fit_err = stdw(d_obs[ind]-d_model[ind],w_d[ind]) ./ (sum(w_d[ind].*dyr[ind])./sum(w_d[ind])); # asg replaced call to wmean [!!! FOUND AND FIXED ERROR !!!!!!]

    t_fit  = Dates.DateTime.(round.(Int,yr),7,1)
    v_fit = p[2+1:end];

elseif model == "interannual"
    amp_fit = missing; # amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).
    phase_fit = missing; # phase converted such that it reflects the day when value is maximized

    # A_err is the *velocity* (not displacement) error, which is the displacement error divided by the weighted mean dt:
    amp_fit_err = missing
    ind = valid;
    amp_fit_err = missing;

    t_fit  = Dates.DateTime.(round.(Int,yr),7,1)
    v_fit = p[1:end];
end

# Number of equivalent image pairs per year: (1 image pair equivalent means a full year of data. It takes about 23 16-day image pairs to make 1 year equivalent image pair.)
fit_count = sum(M[valid,:].>0, dims=1);

v_fit_err =  vec(1 ./ sqrt.(sum(w_v[valid].*M[valid,:], dims=1)));

# as per convention, set missing values as non outliers 
outlier[ismissing.(vec(v))] .= false

return t_fit, v_fit, amp_fit, phase_fit, v_fit_err, amp_fit_err, fit_count, fit_outlier_frac, outlier

end