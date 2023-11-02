push!(LOAD_PATH, "src/")
using Documenter, DocumenterCitations, BayesJ

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

makedocs(
    sitename="BayesJ",
    pages=[
        "Home" => "index.md",
        "Setup" => "setup.md",
        "Usage" => "usage.md",
        "Cluster Models" => "models.md",
        "API" => [
            "Public" => "public.md",
            "Internals" => "private.md"
        ],
        "References" => "references.md"
    ],
    plugins=[bib],
)

deploydocs(
    repo="github.com/b.com/Infinite-Improbability/BayesX-Julia.git",
)