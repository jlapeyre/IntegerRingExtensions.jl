using RingExtensionsCommon
using Documenter
using DocumenterInterLinks: DocumenterInterLinks

DocMeta.setdocmeta!(RingExtensionsCommon, :DocTestSetup, :(using RingExtensionsCommon); recursive=true)

# Example of extref:  [`Statistics.mean`](@extref Julia)

# links = DocumenterInterLinks.InterLinks(
#     "Distributions" => "https://juliastats.org/Distributions.jl/stable/",
#     "StasBase" => "https://juliastats.org/StatsBase.jl/stable/",
#     "Julia" => "https://docs.julialang.org/en/v1/",
# );

makedocs(;
    modules=[RingExtensionsCommon],
    authors="John Lapeyre",
    sitename="RingExtensionsCommon.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:none,
)
