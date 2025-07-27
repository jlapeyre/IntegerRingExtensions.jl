[![Build Status](https://github.com/jlapeyre/IntegerRingExtensions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jlapeyre/IntegerRingExtensions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jlapeyre/IntegerRingExtensions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/IntegerRingExtensions.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

# IntegerRingExtensions

This package provides some extensions of the ring of integers, and tools for working with 2 x 2 matrices over these rings.
Implementations and tools for one-qubit Clifford + T gates are included.

These packages are not registered.

These packages began as tools to verify results from other rotation synthesis software and to verify
calculations appearing in papers.
Routines for doing synthesis are *not* included here. However, these libraries are a suitable base for writing such routines.
Furthermore, this software includes wrapper files [gridsynth.jl](./src/gridsynth.jl) and
[gridsynth_extra.jl](./src/gridsynth_extra.jl) providing a Julia interface to the
original [gridsynth](https://www.mathstat.dal.ca/~selinger/newsynth/) software (written in Haskell).

### Organization

This monorepo contains several installable packages (in [./lib](./lib)).
In addition, there is some [code at the top level](./src) which has not yet been devolved into packages.

These packages are not registered.

The top-level [test suite](./test) runs tests for both the toplevel code as well as the test suits
for each individual package.

* This software is available "as is". I may, or may not, develop it further.
* The level of documentation and testing is uneven.
* Some interfaces are more experimental than others and may be incomplete.

### Related packages

* [Cyclotomics.jl](https://github.com/kalmarek/Cyclotomics.jl)
* [CyclotomicNumbers.jl](https://github.com/jmichel7/CyclotomicNumbers.jl)

These packages implement more general cyclotomic numbers. They are field extensions, rather than
just ring extensions. The implementations are necessarily more complicated and are much less efficient than
those in `IntegerRingExtensions.jl`. In particular, the implementations in these two packages require heap allocation
for creating and using the provided types, whereas instances of rings in `IntegerRingExtensions` are `isbits`
when the chosen parametric integer type is `isbits`.

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlapeyre.github.io/IntegerRingExtensions.jl/dev/) -->


<!--  LocalWords:  Nemo jl ZZRingElem BigInt Matrix2x2 IntegerRingExtensions one-qubit Dev
<!--  LocalWords:  Cyclotomics CyclotomicNumbers cyclotomic ColPrac
 -->

