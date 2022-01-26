# the first time this is run you need to add packages
# using Pkg
# Pkg.add("Zarr")
# Pkg.add("ESDL")
# Pkg.add("ArchGDAL")
# Pkg.add("AWS")
# Pkg.add("YAXArrays")
# Pkg.add("Statistics")
# Pkg.add("Dates")
# Pkg.add("Plots")
# Pkg.add("IterTools")
# Pkg.add("Distributed")
# Pkg.add("JuliennedArrays")
# Pkg.add(PackageSpec(url="https://github.com/esa-esdl/ESDL.jl"))
# Pkg.add("DateFormats")

# load in libraries
using AWS, Zarr, Statistics, Dates, Plots, IterTools, Distributed, JuliennedArrays, ESDL, DateFormats

# set up aws configuration
AWS.global_aws_config(AWSConfig(creds=nothing, region = "us-west-2"))

# for converting numpy date to DateTime (days2sec = 86400)
npdt64_to_dt(time) =  Second.(round.(time.*86400)) .+ DateTime(1970) 

# path to datacube
path = Array{String,1}()
#push!(path,"http://its-live-data.s3.amazonaws.com/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")
#push!(path,"s3://its-live-data/datacubes/v02/N60W040/ITS_LIVE_vel_EPSG3413_G0120_X-150000_Y-2250000.zarr")
push!(path,"s3://its-live-data/datacubes/v02/N60W130/ITS_LIVE_vel_EPSG3413_G0120_X-3250000_Y250000.zarr")
function threaded_read(xin)
    xout = similar(xin)
    Threads.@threads for i in map(i->i.indices,Zarr.DiskArrays.eachchunk(xin))
        xout[i...] = xin[i...]
    end
    xout
end


i = 1;

# map into datacube
z = zopen(path[i], consolidated=true)

# map to specific variables
vx = z["vx"]
vy = z["vy"]

# read 1D variables into memory
dt = threaded_read(z["date_dt"]) 
stable_shift = threaded_read(z["vx_stable_shift"]) 
vx_error = threaded_read(z["vx_error"]) 
vy_error = threaded_read(z["vy_error"]) 
mid_date = threaded_read(z["mid_date"])

# convert data from numpy64 to Julia DateTime
mid_date = npdt64_to_dt(mid_data)

# add systimatic error based on level of co-registration
vy_error[stable_shift.==0] .= vy_error[stable_shift.==0] .+ 100
vx_error[stable_shift.==2] .= vx_error[stable_shift.==2] .+ 20
vy_error[stable_shift.==2] .= vy_error[stable_shift.==2] .+ 20
vx_error[stable_shift.==1] .= vx_error[stable_shift.==1] .+ 5
vy_error[stable_shift.==1] .= vy_error[stable_shift.==1] .+ 5

# Loop for each vertical column of chubked arrays
cs = Zarr.DiskArrays.eachchunk(vx).chunksize[1:2]
sz = size(vx,[1,2])
xout = zeros(sz)
stpi = cld(cs[1],1)
stpj = cld(cs[2],1)
nanmed(x) = median(skipmissing(x))

# read in image-pair time seperation 
dt_edge = [0 32 64 128 256 inf]; # define edges of dt bins

function dtmaxfind(x::Base.ReshapedArray{Union{Missing, Int16}, 3, Array{Union{Missing, Int16}, 3}, Tuple{}}, dt::Vector{Float64}, dt_max::Matrix{Int64})

    # define a median absolut difference function
    function medmad(x::Array{Union{Missing, Int16}}) 
        if all(ismissing, x)
            medx = 0
            madx = 0
        else
            medx = round(Int,median(skipmissing(x)));
            madx = round(Int,median(skipmissing(abs.(x .- medx))));
        end
        return [medx, madx]
    end

    # find last valid bin or return zero
    function findlastorzero(x::BitVector)
        if any(x)
            x = findlast(x)
        else
            x = 0
        end
        return x
    end

    dt = round.(Int, dt)

    # sort index
    dt_sortperm = sortperm(dt);
    dt = dt[dt_sortperm]

    # find the last sorted value that is < 
    bin_ind = searchsortedlast.(Ref(dt), dt_max) 
    bin_ind = [1 bin_ind]

    x = x[:,:,dt_sortperm]

    binMad = zeros(Int16, size(x,1), size(x,2), length(dt_max))
    binMed = zeros(Int16, size(x,1), size(x,2), length(dt_max))

    size(binMad)
    for i in range(1, length(bin_ind)-1)
            foo = mapslices(medmad,x[:,:,bin_ind[i]:bin_ind[i+1]],dims=3);
            binMed[:,:,i] = foo[:,:,1]
            binMad[:,:,i] = foo[:,:,2]
    end

    # dtbin_mad_thresh: used to determine in dt means are significantly different
    dtbin_mad_thresh = 2; 

    # check if populations overlap (use first, smallest dt, bin as reference)
    minBound = binMed - binMad * dtbin_mad_thresh * 1.4826;
    maxBound = binMed + binMad * dtbin_mad_thresh * 1.4826;

    exclude = (minBound .> maxBound[:,:,1]) .| (maxBound .< minBound[1]);

    maxdt = mapslices(findlastorzero, exclude, dims = 3)
    noMax = maxdt .== 0
    maxdt[.~noMax] = dt_max[maxdt[.~noMax]]
    maxdt[noMax] .= (2^15 - 1)

    return maxdt
end

function [A,ph,A_err,t_int,v_int,v_int_err,N_int,outlier_frac] = itslive_lsqfit_annual(t1,t2,v,v_err,mad_thresh,mad_filt_iterations)
   

    ## This code is adapted from Chad A. Greene's Matlab code itslive_lsqfit_annual.m
    # Chad A. Greene, Jan 2020.
    #
    
    # Ensure we're starting with finite data:
    
    # this function fit's 
    # [A,ph,A_err,t_int,v_int,v_int_err,N_int,outlier_frac] 
    
    valid = .~ismissing.(v) 

    t1 = t1[valid];
    t2 = t2[valid];
    v = v[valid];
    v_err = v_err[valid];
    
    
    # Bookkeeping:
    
    # Convert datenums to decimal years:
    yr1 = YearDecimal.(t1)
    yr2 = YearDecimal.(t2)

    # dt in years:
    dyr = yr2 .- yr1
    
    # weights for velocities:
    w_v = 1./v_err.^2;
    
    % Weights (correspond to displacement error, not velocity error):
    w_d = 1./(v_err.*dyr); % Not squared because the p= line below would then have to include sqrt(w) on both accounts
    
    d_obs = v.*dyr; % observed displacement in meters
    
    %% Initial outlier detection and removal:
    % % % I've tried this entire function with or without this section, and it
    % % % doesn't seem to make much difference. It's actually fine to do
    % % % the whole analysis without detrending at all, but I'm leaving this here
    % % % for now anyways.
    % %
    % % % Define outliers as anything that's more than 10 standard deviations away from detrended vals.
    % % outliers = abs(v) > 10*std(v);
    % %
    % % % Remove them!
    % % yr = yr(~outliers,:);
    % % dyr = dyr(~outliers);
    % % d_obs = d_obs(~outliers);
    % % w_d = w_d(~outliers);
    % % w_v = w_v(~outliers);
    
    %% Make matrix of percentages of years corresponding to each displacement measurement
    
    y1 = floor(min(yr(:,1))):floor(max(yr(:,2)));
    
    M = zeros(length(dyr),length(y1));
    
    % Loop through each year:
    for k = 1:length(y1)
        
        % Set all measurements that begin before the first day of the year and end after the last
        % day of the year to 1:
        ind = yr(:,1)<=y1(k) & yr(:,2)>=(y1(k)+1);
        M(ind,k) = 1;
        
        % Within year:
        ind = yr(:,1)>=y1(k) & yr(:,2)<(y1(k)+1);
        M(ind,k) = dyr(ind);
        
        % Started before the beginning of the year and ends during the year:
        ind = yr(:,1)<y1(k) & yr(:,2)>=y1(k) & yr(:,2)<(y1(k)+1);
        M(ind,k) = yr(ind,2)-y1(k);
        
        % Started during the year and ends the next year:
        ind = yr(:,1)>=y1(k) & yr(:,1)<(y1(k)+1) & yr(:,2)>=(y1(k)+1);
        M(ind,k) = (y1(k)+1)-yr(ind,1);
    end
    
    hasdata = sum(M)>0;
    y1 = y1(hasdata);
    M = M(:,hasdata);
    
    %% Iterative mad filter
    totalnum = length(yr);
    for i = 1:mad_filt_iterations
        % Displacement Vandermonde matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
        D = [(cos(2*pi*yr(:,1)) - cos(2*pi*yr(:,2)))./(2*pi) (sin(2*pi*yr(:,2)) - sin(2*pi*yr(:,1)))./(2*pi) M ones(size(dyr))];
        
        % Solve for coefficients of each column in the Vandermonde:
        p = (w_d.*D) \ (w_d.*d_obs);
        
        %% Find and remove outliers
        
        d_model = sum(bsxfun(@times,D,p'),2); % modeled displacements (m)
        
        d_resid = abs(d_obs - d_model)./dyr; % devide by dt to avoid penalizing long dt [asg]
        
        d_sigma = median(d_resid)*1.4826; % robust standard deviation of errors, using median absolute deviation
        
        outliers = d_resid > (mad_thresh*d_sigma);
        
        % Remove outliers:
        yr = yr(~outliers,:);
        dyr = dyr(~outliers);
        d_obs = d_obs(~outliers);
        w_d = w_d(~outliers);
        w_v = w_v(~outliers);
        M = M(~outliers,:);
        
        % Remove no-data columns from M:
        hasdata = sum(M)>1;
        y1 = y1(hasdata);
        M = M(:,hasdata);
        
    end
    outlier_frac = length(yr)./totalnum;
    
    %% Second iteration
    
    % Displacement Vandermonde matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
    D = [(cos(2*pi*yr(:,1)) - cos(2*pi*yr(:,2)))./(2*pi).*(M>0) (sin(2*pi*yr(:,2)) - sin(2*pi*yr(:,1)))./(2*pi).*(M>0) M];
    
    % Solve for coefficients of each column in the Vandermonde:
    p = (w_d.*D) \ (w_d.*d_obs);
    
    %% Postprocess
    
    % Convert coefficients to amplitude and phase of a single sinusoid:
    Nyrs = length(y1);
    A = hypot(p(1:Nyrs),p(Nyrs+1:2*Nyrs)); % amplitude of sinusoid from trig identity a*sin(t) + b*cos(t) = d*sin(t+phi), where d=hypot(a,b) and phi=atan2(b,a).
    
    ph_rad = atan2(p(Nyrs+1:2*Nyrs),p(1:Nyrs)); % phase in radians
    ph = 365.25*(mod(0.25 - ph_rad/(2*pi),1)); % phase converted such that it reflects the day when value is maximized
    
    if nargout>2
        % Goodness of fit:
        d_model = sum(bsxfun(@times,D,p'),2); % modeled displacements (m)
        
        % A_err is the *velocity* (not displacement) error, which is the displacement error divided by the weighted mean dt:
        A_err = NaN(size(A));
        for k = 1:Nyrs
            ind = M(:,k)>0;
            A_err(k) = std(d_obs(ind)-d_model(ind),w_d(ind))./(sum(w_d(ind).*dyr(ind))./sum(dyr(ind))); % asg replaced call to wmean
        end
        
        if nargout>3
            t_int = datenum(y1+0.5,0,0)';
            
            v_int = p(2*Nyrs+1:end);
            % % v_int = p(2*Nyrs+1:end) + polyval(pv,t_int,S,mu);
            
            % Number of equivalent image pairs per year: (1 image pair equivalent means a full year of data. It takes about 23 16-day image pairs to make 1 year equivalent image pair.)
            N_int = sum(M>0);
            
            v_int_err =  1./sqrt(sum(w_v.*M))';
            
        end
    end
    
    end

# define time edges of composites
startDate = 2013;
endDate = 2021;
yrs = startDate:endDate;

# initialize variables
v = Array{Union{Nothing, Float32}}(nothing, sz[1], sz[2], length(yrs))
vx = copy(v)
vy = copy(v)

v_err = copy(v)
vx_err = copy(v)
vy_err = copy(v)

v_cnt = copy(v)
vx_cnt = copy(v)
vy_cnt = copy(v)

v_amp = copy(v)
vx_amp = copy(v)
vy_amp = copy(v)

v_phase = copy(v)
vx_phase = copy(v)
vy_phase = copy(v)

v_sigma = copy(v)
vx_sigma = copy(v)
vy_sigma = copy(v)

outlier_frac = Array{Union{Nothing, Float32}}(nothing, sz[1], sz[2])
maxdt = copy(outlier_frac)

_, fileName = splitdir(path[i])

maxdt = zeros(Int16,sz)
# Threads.@threads for i in range(1, step=stpi, length=cld(sz[1],stpj))
for i in range(1, step=stpi, length=cld(sz[1],stpj))  
    j = for j in range(1, step=stpi, length=cld(sz[2],stpj))
        println([i,j])
 
        # Threaded read hangs - I think this has to do with the choice zarr chunk compression  
        # v2 = view(v, i:min(sz[1], i*stpi), j:min(sz[2], j*stpj), :)
       
        # threaded load into memory
        # @time v2 = threaded_read(v2)

        vxin = view(vx, i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :)
        @time vxin = vxin[:,:,:]

        # vxin = vx[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :]
    
        @time maxdt_vx = dtmaxfind(vxin, dt, dt_max)

        vyin = view(vy, i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :)
        vyin = vyin[:,:,:]

        # vyin = vy[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1), :]
        maxdt_vy = dtmaxfind(vyin, dt, dt_max)

        foo = cat(maxdt_vx, maxdt_vy,dims = 3)
        maxdt[i:min(sz[1], i+stpi-1), j:min(sz[2], j+stpj-1)] = minimum(foo,dims=3)



        return
    end
end

gr()
heatmap(1:size(maxdt,1),
    1:size(maxdt,2), maxdt,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="Median Velocity")