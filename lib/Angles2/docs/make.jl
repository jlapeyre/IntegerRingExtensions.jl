using Angles2
using Documenter

DocMeta.setdocmeta!(Angles2, :DocTestSetup, :(using Angles2); recursive=true)

makedocs(;
    modules=[Angles2],
    authors="John Lapeyre",
    sitename="Angles2.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs = :none,
)
