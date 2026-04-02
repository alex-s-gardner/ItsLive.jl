# Feature request: native URL/HTTP support and better error guidance in `GeoDataFrames.read`

## Summary

`GeoDataFrames.read` I wonder if GeoDataFrames could better support remote files... it's quite possible I'm using GeoDataFrames incorrectly

## Background

Below is what I think a typical (non-expert) user would attempt to do:

```julia
using GeoDataFrames, HTTP, GeoJSON

path = "https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json"

# Attempt 1: direct URL string
GeoDataFrames.read(path)
```

<details>
<summary>Cryptic error output</summary>

```
┌ Warning: Failed to parse GeoJSON as 2D, trying 3D. Set `ndim` to 3 to avoid this warning.
└ @ GeoJSON ~/.julia/packages/GeoJSON/73eSw/src/io.jl:21
ERROR: ArgumentError: invalid JSON at byte position 1 while parsing type GeoJSON.GeoJSONT{3, Float32}: ExpectedOpeningObjectChar
https://its-live-data.s3-u

Stacktrace:
  [1] invalid(error::JSON3.Error, buf::Base.CodeUnits{UInt8, String}, pos::Int64, T::Type)
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/JSON3.jl:30
  [2] parse(str::String, ::Type{GeoJSON.GeoJSONWrapper{…} where X<:GeoJSON.GeoJSONT{…}}; jsonlines::Bool, kw::@Kwargs{})
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/structs.jl:32
  [3] parse
    @ ~/.julia/packages/JSON3/rT1w2/src/structs.jl:32 [inlined]
  [4] #read#18
    @ ~/.julia/packages/JSON3/rT1w2/src/structs.jl:46 [inlined]
  [5] read(str::String, ::Type{GeoJSON.GeoJSONWrapper{3, Float32, X} where X<:GeoJSON.GeoJSONT{3, Float32}})
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/structs.jl:45
  [6] read(io::String; lazyfc::Bool, ndim::Int64, numbertype::Type)
    @ GeoJSON ~/.julia/packages/GeoJSON/73eSw/src/io.jl:22
  [7] read
    @ ~/.julia/packages/GeoJSON/73eSw/src/io.jl:13 [inlined]
  [8] read(::GeoDataFrames.GeoJSONDriver, fname::String; kwargs::@Kwargs{})
    @ GeoDataFramesGeoJSONExt ~/.julia/packages/GeoDataFrames/YTFj5/ext/GeoDataFramesGeoJSONExt.jl:17
  [9] read
    @ ~/.julia/packages/GeoDataFrames/YTFj5/ext/GeoDataFramesGeoJSONExt.jl:13 [inlined]
 [10] read(fn::String; kwargs::@Kwargs{})
    @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:53
 [11] read(fn::String)
    @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:50
 [12] top-level scope
    @ REPL[1]:1

caused by: ArgumentError: invalid JSON at byte position 1 while parsing type GeoJSON.GeoJSONT{2, Float32}: ExpectedOpeningObjectChar
https://its-live-data.s3-u

Stacktrace:
  [1] invalid(error::JSON3.Error, buf::Base.CodeUnits{UInt8, String}, pos::Int64, T::Type)
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/JSON3.jl:30
  [2] parse(str::String, ::Type{GeoJSON.GeoJSONWrapper{…} where X<:GeoJSON.GeoJSONT{…}}; jsonlines::Bool, kw::@Kwargs{})
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/structs.jl:32
  [3] parse
    @ ~/.julia/packages/JSON3/rT1w2/src/structs.jl:32 [inlined]
  [4] #read#18
    @ ~/.julia/packages/JSON3/rT1w2/src/structs.jl:46 [inlined]
  [5] read(str::String, ::Type{GeoJSON.GeoJSONWrapper{2, Float32, X} where X<:GeoJSON.GeoJSONT{2, Float32}})
    @ JSON3 ~/.julia/packages/JSON3/rT1w2/src/structs.jl:45
  [6] read(io::String; lazyfc::Bool, ndim::Int64, numbertype::Type)
    @ GeoJSON ~/.julia/packages/GeoJSON/73eSw/src/io.jl:18
  [7] read
    @ ~/.julia/packages/GeoJSON/73eSw/src/io.jl:13 [inlined]
  [8] read(::GeoDataFrames.GeoJSONDriver, fname::String; kwargs::@Kwargs{})
    @ GeoDataFramesGeoJSONExt ~/.julia/packages/GeoDataFrames/YTFj5/ext/GeoDataFramesGeoJSONExt.jl:17
  [9] read
    @ ~/.julia/packages/GeoDataFrames/YTFj5/ext/GeoDataFramesGeoJSONExt.jl:13 [inlined]
 [10] read(fn::String; kwargs::@Kwargs{})
    @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:53
 [11] read(fn::String)
    @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:50
 [12] top-level scope
    @ REPL[1]:1
Some type information was truncated. Use `show(err)` to see complete types.
```

</details>

```julia
# Attempt 2: explicit GDAL driver with URL
GeoDataFrames.read(GeoDataFrames.ArchGDALDriver(), path)
```

<details>
<summary>Cryptic error output</summary>

```
GDALError (CE_Failure, code 1):
	SSL certificate problem: unable to get local issuer certificate

Stacktrace:
  [1] maybe_throw
    @ ~/.julia/packages/GDAL/3ZrbT/src/error.jl:42 [inlined]
  [2] aftercare
    @ ~/.julia/packages/GDAL/3ZrbT/src/error.jl:59 [inlined]
  [3] gdalopenex(pszFilename::String, nOpenFlags::Int64, ...)
    @ GDAL ~/.julia/packages/GDAL/3ZrbT/src/libgdal.jl:7540
  [4] read(driver::GeoDataFrames.ArchGDALDriver, fn::String)
    @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:84
```

</details>

```julia
# Attempt 3: /vsicurl/ virtual filesystem prefix (standard GDAL convention for remote files)
GeoDataFrames.read(GeoDataFrames.ArchGDALDriver(), "/vsicurl/" * path)
```

<details>
<summary>Cryptic error output</summary>

```
ERROR: Unable to open /vsicurl/https://its-live-data.s3-us-west-2.amazonaws.com/datacubes/catalog_v02.json.
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:44
 [2] (::GeoDataFrames.var"#6#7"{Nothing, GeoDataFrames.ArchGDALDriver, String})(ds::ArchGDAL.Dataset)
   @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:87
 [3] read(f::GeoDataFrames.var"#6#7"{Nothing, GeoDataFrames.ArchGDALDriver, String}, args::String; kwargs::@Kwargs{})
   @ ArchGDAL ~/.julia/packages/ArchGDAL/7ZWOl/src/context.jl:268
 [4] read(driver::GeoDataFrames.ArchGDALDriver, fn::String)
   @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:84
 [5] top-level scope
   @ REPL[3]:1
```

</details>

```julia
# Attempt 4: pass HTTP.Response object directly
gj = HTTP.get(path)
GeoDataFrames.read(gj)
```

<details>
<summary>Cryptic error output</summary>

```
ERROR: MethodError: no method matching splitext(::HTTP.Messages.Response)
The function `splitext` exists, but no method is defined for this combination of argument types.

Closest candidates are:
  splitext(::String)
   @ Base path.jl:226
  splitext(::AbstractString)
   @ Base path.jl:614

Stacktrace:
 [1] read(fn::HTTP.Messages.Response; kwargs::@Kwargs{})
   @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:51
 [2] read(fn::HTTP.Messages.Response)
   @ GeoDataFrames ~/.julia/packages/GeoDataFrames/YTFj5/src/io.jl:50
 [3] top-level scope
   @ REPL[5]:1
```

</details>

```julia
# Attempt 5: pass the HTTP response body as a raw GeoJSON string — this works!
gj = String(HTTP.get(path).body)
GeoDataFrames.read(gj)
# => Returns a 3160×10 DataFrame
```

## Environment
- GeoDataFrames: v0.4.2
- GeoJSON: v0.8.4