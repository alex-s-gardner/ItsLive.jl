"""
    annual_matrix(t1,t2)

create an annual weighting matrix for discrete intervals

# Example
```julia
julia> M, tM = annual_matrix(t1,t2)
```

# Arguments
   - `t1::Vector{DateTime}`: start DateTime of descrete interval 
   - `t2::Vector{DateTime}`: end DateTime of descrete interval 

# Author
Alex Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 23, 2022

Chad A. Greene [inspiration Matlab code]
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
January 1, 2022
"""
function annual_matrix(t1::Vector{DateTime}, t2::Vector{DateTime})
    # Convert datenums to decimal years:
    yr1 = ItsLive.decimalyear(t1)
    yr2 = ItsLive.decimalyear(t2)
    
    # make annual weighting matrix for discrete intervals
    # unique years
    yr = floor(minimum(yr1)):floor(maximum(yr2));

    # dt in years:
    dyr = yr2 .- yr1

    M = zeros(length(yr1),length(yr));

    ## Make matrix of percentages of years corresponding to each displacement measurement

    # Loop through each year:
    for k = 1:length(yr)
        
        # Set all measurements that begin before the first day of the year and end after the last
        # day of the year to 1:
        ind = (yr1.<=yr[k]) .& (yr2.>=(yr[k]+1));
        M[ind,k] .= 1;
        
        # Within year:
        ind = (yr1.>=yr[k]) .& (yr2.<(yr[k]+1));
        M[ind,k] = dyr[ind];
        
        # Started before the beginning of the year and ends during the year:
        ind = (yr1.<yr[k]) .& (yr2.>=yr[k]) .& (yr2.<(yr[k].+1));
        M[ind,k] = yr2[ind] .- yr[k];
        
        # Started during the year and ends the next year:
        ind = (yr1.>=yr[k]) .& (yr1.<(yr[k].+1)) .& (yr2.>=(yr[k]+1));
        M[ind,k] = (yr[k]+1) .- yr1[ind];
    end
    tM  = Dates.DateTime.(round.(Int,yr),7,1)

    return M, tM
end