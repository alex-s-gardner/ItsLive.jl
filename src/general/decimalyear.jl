"""
    decyear = decimalyear(datetime)

return the decimal year (`decyear`) of a DateTime type vector

# Example
````julia
    julia> decyear = ItsLive.decimalyear(Dates.DateTime(1970,1,1))
    1970.0
````

# Arguments
   - `datetime::Union{DateTime, Vector{DateTime}}`: date and time of type DateTime 

# Author
Alex S. Gardner, JPL, Caltech.
"""
function decimalyear(datetime::Union{DateTime, Vector{DateTime}})
     
     # function for converting YearDecimal
     function YearDecimal2Float64(x)
        decyr = Vector{Union{Missing, Float64}}(undef,length(x))
        for i = 1:length(x) 
            decyr[i] = x[i].value 
        end
        return decyr
    end

    # Convert datenums to decimal years:
    yr1 = DateFormats.YearDecimal.(datetime)
    decimalyr = YearDecimal2Float64(yr1)
    return decimalyr
end

function decimalyear(datetime::DateTime)

    decimalyr = DateFormats.YearDecimal.(datetime)
    decimalyr = decimalyr.value 

   return decimalyr
end
