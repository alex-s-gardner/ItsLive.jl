using Documenter
using ItsLive

makedocs(
    sitename = "ItsLive",
    format = Documenter.HTML(),
    modules = [ItsLive]
)

deploydocs(
    repo = "github.com/alex-s-gardner/ItsLive.jl.git",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
