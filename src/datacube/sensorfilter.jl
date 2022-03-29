"""
    sensorfilter(vx, vy, mid_date, dt, sensor; binedges, mincount, dtmax, id_refsensor, plotflag)

    identify sensors that produce velocities that are sginificanlty slower

# Example no inputs
```julia
julia> sensor_exclude = sensorfilter(C[1,"vx"], C[1,"vy"], C[1,"mid_date"], C[1,"date_dt"], C[1,"satellite_img1"])
```

# Arguments
    - `vx::Vector{Any}`: x component velocity
    - `vy::Vector{Any}`: y component velocity
    - `mid_date::Vector{DateTime}`: center date of image-pair [] 
    - `dt::Vector{Any}`: time seperation between image pairs [days]
    - `sensor::Vector{Any}`: list of image sensors for image-pairs... used to group results

# Keyword Arguments
    - `binedges::Union{Vector{Number}, StepRangeLen{Number, Base.TwicePrecision{Number}, Base.TwicePrecision{Number}, Number}}`: DateTime bin edges
    - `mincount::Number`: minimum count for a bin to be considered valid
    - `dtmax::Number`: maximum dt in days of admissible data for comparison
    - `id_refsensor::Number`: id of sensor to which all other sensores are compared
    - `sescale::Number`: standard error (se) scale factor for determination of significance
    - `plotflag::Bool`: flag for turning on and off plotting

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
March 28, 2022
"""

function sensorfilter(vx::Vector{Union{Missing, Float64}}, vy::Vector{Union{Missing, Float64}},
    mid_date::Vector{DateTime}, dt::Vector{Float64}, sensor::Vector{Any};
    binedges = 2015.5:(1/5):2022.5,
    mincount::Number = 3,
    dtmax::Number = 24,
    id_refsensor::Number = 1,
    sescale::Number = 1,
    plotflag::Bool = true)

    # expand bin edges if needed
    if (binedges isa StepRangeLen)
        binedges = collect(binedges)
    end

    # trim data to redude computations
    dtind = dt .<= dtmax;
    valid = .~ismissing.(vx)
    ind = dtind .& valid

    vx = vx[ind];
    vy = vy[ind];
    sensor = sensor[ind]
    dt = dt[ind]
    mid_date = mid_date[ind]

    # determine sensor group ids
    id, sensorgroups = ItsLive.sensorgroup(sensor)
    numsg = length(sensorgroups)

    # convert date to decimal year
    decyear = ItsLive.decimalyear(mid_date)

    # initialize veriables
    vbin = fill(NaN, (numsg,length(binedges)-1))
    vstdbin = copy(vbin)
    vcountbin = zeros(Int, numsg, length(binedges)-1)
    bincenters = []

    # loop for each sensor group
    for sg = 1:numsg
        ind = id .== sg;
        vx0 = mean(vx[ind])
        vy0 = mean(vy[ind])
        v0 = sqrt.(vx0.^2 .+ vy0.^2);
        uv = vcat(vx0/v0, vy0/v0)
        vp = hcat(vx[ind],vy[ind]) * uv # flow acceleration in direction of unit flow vector
        
        vbin[sg,:], vstdbin[sg,:], vcountbin[sg,:], bincenters = ItsLive.binstats(decyear[ind], vp; binedges)
    end

    # remove points with < mincount
    vbin[vcountbin .< mincount] .= NaN
    
    # plot binned values
    if plotflag
        plot()
        for sg = 1:numsg
            if sum(.~isnan.(vbin[sg,:]))>3
                p = plot!(bincenters, vbin[sg,:]; markershape=:circle, markersize = 4, label = sensorgroups[sg]["name"])
                display(p)    
            end 
        end
    end

    # check if Setinel-2 mean is different from Landsat-8
    # calculate mean (m) and standard error (s)
    m = fill(NaN, numsg)
    s = fill(NaN, numsg)
    for sg = 1:numsg
        covalid = .~isnan.(vbin[id_refsensor,:]) .& .~isnan.(vbin[sg,:])
        if sum(covalid) > 3
            delta =   vbin[sg,covalid] - vbin[id_refsensor,covalid];
            m[sg,1] = mean(delta)
            s[sg,1] = std(delta)/sqrt((sum(covalid)-1))
        end
    end

    # check if the mean difference + se < zero, 
    # if < zero then id_refsensor has a sginificanlty faster mean
    disagree = (m .+ (s.*sescale)) .< 0

    # retrieve sensor id
    sensor_exclude = collect(1:numsg)[disagree]

    return sensor_exclude
end