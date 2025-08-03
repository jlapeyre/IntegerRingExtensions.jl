# IntegerExtensions

This package provides some extensions of the ring of integers and tools for working with 2 x 2 matrices over these rings.
Implementations and tools for one-qubit Clifford + T gates are included.

The main goal of this package is to support research in unitary synthesis.

### Related packages

* [Cyclotomics.jl](https://github.com/kalmarek/Cyclotomics.jl)
* [CyclotomicNumbers.jl](https://github.com/jmichel7/CyclotomicNumbers.jl)

These packages implement more general cyclotomic numbers. They are field extensions, rather than
just ring extensions. The implementations are necessarily more complicated and are much less efficient than
those in `IntegerExtensions.jl`. In particular, the implementations in these packages require heap allocation
for creating and using the provided types.

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlapeyre.github.io/IntegerExtensions.jl/dev/) -->


<!-- [![Build Status](https://github.com/jlapeyre/IntegerExtensions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jlapeyre/IntegerExtensions.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

<!-- [![Coverage](https://codecov.io/gh/jlapeyre/IntegerExtensions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/IntegerExtensions.jl) -->

<!-- [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) -->


<!--  LocalWords:  Nemo jl ZZRingElem BigInt Matrix2x2 IntegerExtensions one-qubit Dev
<!--  LocalWords:  Cyclotomics CyclotomicNumbers cyclotomic ColPrac
 -->

