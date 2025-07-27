using SingletonNumbers
using Documenter

DocMeta.setdocmeta!(SingletonNumbers, :DocTestSetup, :(using SingletonNumbers); recursive=true)

makedocs(;
    modules=[SingletonNumbers],
    authors="John Lapeyre",
    sitename="SingletonNumbers.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
         checkdocs = :none,
         doctest = true,
)
