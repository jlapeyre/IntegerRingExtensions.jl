using IntegerExtensions
using Documenter

DocMeta.setdocmeta!(IntegerExtensions, :DocTestSetup, :(using IntegerExtensions); recursive=true)

makedocs(;
    modules=[IntegerExtensions],
    authors="John Lapeyre",
    sitename="IntegerExtensions.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
