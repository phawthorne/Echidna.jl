using Documenter, Echidna
push!(LOAD_PATH,"../src/")
makedocs(
    format = Documenter.HTML(),
    modules = [Echidna],
    sitename = "Echidna.jl",
    pages = [
        "Home" => "index.md",
    ]
)
