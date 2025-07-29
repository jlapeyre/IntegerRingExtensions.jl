# IntegerExtensions

This package provides extensions of the ring of integers and tools for working with 2 x 2 matrices over these rings.
Implementations and tools for one-qubit Clifford + T gates are included.

This package exists to support research in unitary synthesis.

### Related packages

* [Cyclotomics.jl](https://github.com/kalmarek/Cyclotomics.jl)
* [CyclotomicNumbers.jl](https://github.com/jmichel7/CyclotomicNumbers.jl)

These packages implement more general cyclotomic numbers. They are field extensions, rather than
jus ring extensions. The implementations are much less efficient than
those in `IntegerExtensions.jl`. In particular, the implementations in these packages require heap alloction
for creating and using the provided types.

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlapeyre.github.io/IntegerExtensions.jl/stable/) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlapeyre.github.io/IntegerExtensions.jl/dev/) -->


<!-- [![Build Status](https://github.com/jlapeyre/IntegerExtensions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jlapeyre/IntegerExtensions.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

<!-- [![Coverage](https://codecov.io/gh/jlapeyre/IntegerExtensions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/IntegerExtensions.jl) -->

<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->

