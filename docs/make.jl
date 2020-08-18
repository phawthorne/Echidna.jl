# using Documenter, Echidna
# push!(LOAD_PATH,"../src/")
# makedocs(
#     format = Documenter.HTML(),
#     modules = [Echidna],
#     sitename = "Echidna.jl",
#     pages = [
#         "Home" => "index.md",
#     ]
# )

import Pkg
Pkg.activate("..")

using Documenter, Echidna

DocMeta.setdocmeta!(Echidna, :DocTestSetup,
                    :(using Echidna); recursive=true)

makedocs(
    modules = [Echidna],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    clean = true,
    sitename = "Echidna.jl",
    authors = "Peter Hawthorne",
    checkdocs = :exports,
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/phawthorne/Echidna.git"
)