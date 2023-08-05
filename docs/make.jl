push!(LOAD_PATH, "src/")
using Documenter, BayesJ

makedocs(
    sitename="BayesJ",
    pages=[
        "Home" => "index.md",
        "Reference" => [
            "Public" => "public.md",
            "Internals" => "private.md"
        ]
    ]
)

deploydocs(deploydocs(
    repo="github.com/b.com/Infinite-Improbability/BayesX-Julia.git",
))