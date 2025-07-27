using CyclotomicRings
using Documenter

DocMeta.setdocmeta!(CyclotomicRings, :DocTestSetup, :(using CyclotomicRings); recursive=true)

makedocs(;
    modules=[CyclotomicRings],
    authors="John Lapeyre",
    sitename="CyclotomicRings.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
