using QuadraticRings
using Documenter

DocMeta.setdocmeta!(QuadraticRings, :DocTestSetup, :(using QuadraticRings); recursive=true)

makedocs(;
    modules=[QuadraticRings],
    authors="John Lapeyre",
    sitename="QuadraticRings.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
