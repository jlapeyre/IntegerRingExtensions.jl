using RingExtensionsUtils
using Documenter

DocMeta.setdocmeta!(RingExtensionsUtils, :DocTestSetup, :(using RingExtensionsUtils); recursive=true)

makedocs(;
    modules=[RingExtensionsUtils],
    authors="John Lapeyre",
    sitename="RingExtensionsUtils.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
