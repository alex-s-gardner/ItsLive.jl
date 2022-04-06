using Documenter
using ItsLive

# DocMeta.setdocmeta!(ItsLive, :DocTestSetup, :(using ItsLive); recursive=true)

makedocs(
    modules=[ItsLive],
    authors="Alex S. Gardner, JPL, Caltech.",
    repo="https://github.com/JuliaClimate/STAC.jl/blob/{commit}{path}#{line}",
    sitename = "ItsLive.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alex.s.gardner.github.io/ItsLive.jl/dev/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],  
)

deploydocs(
    repo = "github.com/alex-s-gardner/ItsLive.jl",
    devbranch="main",
)
