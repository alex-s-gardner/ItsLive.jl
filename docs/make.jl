using Documenter
using ItsLive

makedocs(
    modules=[ItsLive],
    authors="Alex S. Gardner, JPL, Caltech",
    format = Documenter.HTML(
        sidebar_sitename = false,
    )
)

deploydocs(
    repo = "github.com/alex-s-gardner/ItsLive.jl.git",
)