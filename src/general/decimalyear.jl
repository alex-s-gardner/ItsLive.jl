"""
    decimalyear(dt)

returns the decimal year of a DateTime type vector

using Dates

# Example
```julia
julia> decimalyear(Dates.DateTime(1970,5,1))
```

# Arguments
   - `dt::Union{DateTime, Vector{DateTime}}`: date and time of DateTime type

# Author
Alex S. Gardner
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California
February 11, 2022
"""

function decimalyear(dt::Union{DateTime, Vector{DateTime}})
     
     # function for converting YearDecimal
     function YearDecimal2Float64(x)
        decyr = Vector{Union{Missing, Float64}}(undef,length(x))
        for i = 1:length(x) 
            decyr[i] = x[i].value 
        end
        return decyr
    end

    # Convert datenums to decimal years:
    yr1 = DateFormats.YearDecimal.(dt)
    decimalyr = YearDecimal2Float64(yr1)
    return decimalyr
end

function decimalyear(dt::DateTime)

    decimalyr = DateFormats.YearDecimal.(dt)
    decimalyr = decimalyr.value 

   return decimalyr
end
