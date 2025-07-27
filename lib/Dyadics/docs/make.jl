using Dyadics
using Documenter

DocMeta.setdocmeta!(Dyadics, :DocTestSetup, :(using Dyadics); recursive=true)

makedocs(;
    modules=[Dyadics],
    authors="John Lapeyre",
    sitename="Dyadics.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
