using RootOnes
using Documenter

DocMeta.setdocmeta!(RootOnes, :DocTestSetup, :(using RootOnes); recursive=true)

makedocs(;
         modules=[RootOnes],
         authors="John Lapeyre",
         sitename="RootOnes.jl",
         format=Documenter.HTML(;
                                edit_link="main",
                                assets=String[],
                                ),
         pages=[
             "Home" => "index.md",
         ],
         doctest=true,
         checkdocs=:none,
)
