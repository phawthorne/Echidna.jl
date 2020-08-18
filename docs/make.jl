import Pkg

if get(ENV, "CI", nothing) === nothing
    Pkg.activate("..")
end

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

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/phawthorne/Echidna.git"
    )
end