"""
    D, tD, M = design_matrix(t1,t2, model [optional])

create a design matrix (`D`) for discrete interval data

# Example
```julia
julia> D, tD, M = design_matrix(t1,t2, model [optional])
```

# Arguments
   - `t1::Vector{DateTime}`: start DateTime of descrete interval 
   - `t2::Vector{DateTime}`: end DateTime of descrete interval 
   - `model::String = "sinusoidal_interannual"`:

   - `D`: Design matrix
   - `tD`: DateTime centers used in `D`
   - `M` = Annual weighting matrix`:

# Author
Alex S. Gardner and Chad A. Greene, JPL, Caltech.
"""
function design_matrix(t1::Vector{DateTime}, t2::Vector{DateTime}, model::String = "sinusoidal_interannual")
    # Convert datenums to decimal years:
    yr1 = ItsLive.decimalyear(t1)
    yr2 = ItsLive.decimalyear(t2)

    if (model == "sinusoidal_interannual") || (model == "sinusoidal") || (model == "interannual")

    M, tM = ItsLive.annual_matrix(t1,t2)
    yr = ItsLive.decimalyear(tM)

    hasdata = vec(sum(M,dims=1).>0);
    yr = yr[hasdata];
    M = M[:,hasdata];

    if model == "interannual"
        D = M; #hcat(M);
    elseif model == "sinusoidal"
        D = hcat((cos.(2*pi*yr1) - cos.(2*pi*yr2)) ./ (2*pi), (sin.(2*pi*yr2) - sin.(2*pi*yr1))./(2*pi), M);
    else
        # Displacement Design matrix: (these are displacements! not velocities, so this matrix is just the definite integral wrt time of a*sin(2*pi*yr)+b*cos(2*pi*yr)+c.
        N = ((sin.(2*pi*yr2) - sin.(2*pi*yr1))./(2*pi));
        N = mapslices(row->(N).*row,float(M.>0), dims=1);

        P = (cos.(2*pi*yr1) - cos.(2*pi*yr2)) ./ (2*pi);
        P = mapslices(row->(P).*row,float(M.>0), dims=1);

        D = hcat(P, N, M);
    end

    tD  = Dates.DateTime.(round.(Int,yr),7,1)
    return D, tD, M
end
    error("$model is not a design matrix model option")
end