using Documenter, Echidna
push!(LOAD_PATH,"../src/")
makedocs(
    format = :html,
    modules = [Echidna],
    sitename = "Echidna.jl",
    pages = [
        "Home" => "index.md",
    ]
)
