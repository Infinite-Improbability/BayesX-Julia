push!(LOAD_PATH, "src/")
using Documenter, DocumenterCitations, BayesJ

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

makedocs(
    bib,
    sitename="BayesJ",
    pages=[
        "Home" => "index.md",
        "Cluster Models" => "models.md",
        "API" => [
            "Public" => "public.md",
            "Internals" => "private.md"
        ],
        "References" => "references.md"
    ]
)

deploydocs(
    repo="github.com/b.com/Infinite-Improbability/BayesX-Julia.git",
)