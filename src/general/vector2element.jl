"""
    df0 = vector2element(df::DataFrame)

return a DataFrame (`df`) composed of the elements of the input DataFrame (`df`)

using DataFrames

# Example
```julia
julia> df = vect2element(df)
```

# Arguments
   - `df::DataFrame`
"""

function vector2element(df::DataFrame)
    # loop through each column
    var_name = names(df)
    var_nelem = zeros(size(var_name))
    for i in eachindex(var_name)
        var_nelem[i] = sum(length.(df[!, var_name[i]]))
    end

    # what are the number of elements that are greater than one element per row
    nelem = floor.(Int, unique(var_nelem[var_nelem .!= nrow(df)]));
    if length(nelem) > 1
        error("inconsistent number of elements per variable")
    end

    # expand single element rows to match nelem of other rows
    ind = findfirst(var_nelem .== nelem)
    for i in findall(var_nelem .== nrow(df))
        df[!,var_name[i]] = [zeros(size(df[j,var_name[ind]])) .+ df[j,var_name[i]] for j = eachindex(df[!,var_name[i]])]
    end

    df0 = DataFrame()
    for i in eachindex(var_name)
        df0[!,var_name[i]] = reduce(vcat, df[!,var_name[i]])
    end

    return df0
end