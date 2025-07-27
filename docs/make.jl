using IntegerRingExtensions

##
## Each documented module must appear in at least three places
## 1) using MyMod
## 2) in `modules = Module[...` below
## 3) implicitly, in the `pages = [...` The md page will refer to the module
##

#using RingExtensionsUtils

# using Angles2
# using SingletonNumbers
using RingExtensionsCommon
using RootOnes
# using Dyadics
# using QuadraticRings
# using CyclotomicRings

using Documenter

DocMeta.setdocmeta!(IntegerRingExtensions, :DocTestSetup, :(using IntegerRingExtensions); recursive=true)

makedocs(;
    # Despite what the Documenter docs say, documented Modules *must* appear in this list.
    modules = Module[IntegerRingExtensions, RootOnes, RingExtensionsCommon],
    authors="John Lapeyre",
    sitename="IntegerRingExtensions.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Packages" => [
            "packages/root_ones.md",
#            "packages/angles2.md",
            # "packages/cyclotomic_rings.md",
            # "packages/dyadics.md",
            # "packages/quadratic_rings.md",
            "packages/ring_extensions_common.md",
            # "packages/ring_extensions_utils.md",
            # "packages/singleton_numbers.md",
        ]
    ],
    doctest=false,
    checkdocs=:none,
)
