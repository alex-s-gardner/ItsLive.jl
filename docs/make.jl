using Documenter
using ItsLive

# DocMeta.setdocmeta!(ItsLive, :DocTestSetup, :(using ItsLive); recursive=true)

makedocs(
    modules=[ItsLive],
    authors="Alex S. Gardner, JPL, Caltech",
    repo="https://github.com/alex-s-gardner/ItsLive.jl/blob/{commit}{path}#{line}",
    sitename = "ItsLive.jl",
)

deploydocs(
    repo = "github.com/alex-s-gardner/ItsLive.jl",
    devbranch="main",
)
