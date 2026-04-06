


const sensor_groupings = Dict{Int, @NamedTuple{id::Int, name::String, sensors::Vector{String}}}(
    1 => (id=1, name="Sentinel 2",   sensors=["2A", "2B", "2C", "2D"]), # <= id `1` is default reference group for inter-sensor comparison
    2 => (id=2, name="Landsat 8/9",  sensors=["8",  "9" , "10"]),
    3 => (id=3, name="Sentinel 1",   sensors=["1A", "1B", "1C", "1D"]),
    4 => (id=4, name="Landsat 7",    sensors=["7"       ]),
    5 => (id=5, name="Landsat 4/5",  sensors=["4",  "5" ]),
)


"""
    to_datetime(ar; yeardecimal=false)

Convert a DimArray variable with time units to DateTime.

# Arguments
- `ar`: A DimensionalData array with time units metadata
- `yeardecimal::Bool=false`: If true, convert DateTime to year decimal format

# Returns
A DimensionalData array with DateTime values

# Throws
- `InexactError`: If time decoding fails
- `ErrorException`: If date units are not supported
"""
function to_datetime(ar; yeardecimal=false)
    aratts = DimensionalData.metadata(ar)
    if match(r"^(days)|(hours)|(seconds)|(months) since", lowercase(get(aratts, "units", ""))) !== nothing
        timevar = try
            ItsLive.timedecode(ar.data, aratts["units"], lowercase(get(aratts, "calendar", "standard")))
        catch e
            return throw(e)
            if e isa InexactError
                dec = timedecode(ar.data, aratts["units"], lowercase(get(aratts, "calendar", "standard")), prefer_datetime=false)
                if force_datetime
                    round_datetime.(dec)
                else
                    dec
                end
            else
                ar
            end
        end

        if yeardecimal
            timevar = yeardecimal.(timevar)
        end
        timevar = DimensionalData.rebuild(ar, data=timevar)

        return timevar
    else
        error("Date Units not supported")
    end
end


"""
    buffer_unidirectional(line, buffer; dims = 1)

Create a buffered polygon from a line geometry, applying buffer in one direction.

# Arguments
- `line`: A geometry object with a `.geom` attribute containing coordinates
- `buffer`: The buffer distance to apply
- `dims::Int=1`: Dimension to apply buffer (1 for x-direction, 2 for y-direction)

# Returns
A `GI.Polygon` object representing the buffered geometry

# Throws
- `ErrorException`: If dims is not 1 or 2
"""
function buffer_unidirectional(line, buffer; dims = 1)

    cords = GI.coordinates.(line.geom)
    b1 = zero(buffer)
    b2 = zero(buffer)
    if dims == 1
        b1 = buffer
    elseif dims == 2
        b2 = buffer
    else
        error("Invalid dimension. Must be 1 or 2.")
    end

    up = tuple.(getindex.(cords, 1) .+ b1, getindex.(cords, 2) .+ b2)
    dn = reverse(tuple.(getindex.(cords, 1) .- b1, getindex.(cords, 2) .- b2))
    buffered_polygon = GI.Polygon(GI.LineString(GI.Point.(vcat(up, dn))))
    return buffered_polygon
end