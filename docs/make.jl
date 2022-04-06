using Documenter
using ItsLive

makedocs(
    modules=[ItsLive],
    authors="Alex S. Gardner, JPL, Caltech",
    sitename="ItsLive.jl",
    canonical = "https://alex-s-gardner.github.io/ItsLive.jl/stable/",
    format = Documenter.HTML(
        sidebar_sitename = false,
    )
)

deploydocs(
    repo = "github.com/alex-s-gardner/ItsLive.jl.git",
    target = "build",
    branch = "gh-pages",
    forcepush = true,
    dirname = "",
    push_preview = true,
)