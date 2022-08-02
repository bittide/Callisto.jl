
using Documenter
using Callisto




makedocs(sitename="Callisto Documentation", format = Documenter.HTML())

deploydocs(versions=nothing, repo = "github.com/bittide/Callisto.jl.git",)



