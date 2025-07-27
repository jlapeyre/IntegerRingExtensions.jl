module Benchmarking

# This requires the user to import a weak dependency manually.
# BenchmarkTools is slightly heavy compared to the rest of the project.

# See ../ext/BenchmarkToolsExt.jl

export benchmark_compose

"""
    benchmark_compose(gate_str)

Measure and return the execution time per gate in `ns` of `compose(gate_str)`.

To use `benchmark_compose`, you must have `BenchmarkTools` in your load path and `import BenchmarkTools`.

Timing is done via `@belapsed`.

# Examples

Compose the output of `gridsynth pi/8 -f 100 -s -p -d 10`
```julia
julia> using IntegerRingExtensions, BenchmarkTools, Printf;

julia> gates_r_pi_over_8 = "HTHTSHTHTHTSHTHTSHTSHTSHTSHTHTSHTSHTHTSHTSHTHTHTHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTSHTHTSHTHTHTHTHTSHTSHTHTSHTSHTSHTHTHTHTSHTSHTHTSHTSHTHTHTHTSHTSHTSHTSHTSHTSHTSHTSHTHTHTSHTSHTHTHTHTHTHTSHTSHTSHTHTHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTSHTSHTHTSHTHWWWWWWW";

julia> m = compose(gates_r_pi_over_8)
2×2 Matrix2x2{DOmega{Int64}}:
25541547/2²⁵ ω⁰ + 6324583/2²⁵ ω + -4061155/2²⁴ ω² + -4095549/2²⁵ ω³  3294431/2²³ ω⁰ + -4192749/2²⁴ ω + -1318839/2²⁵ ω² + 5125309/2²⁴ ω³
-3294431/2²³ ω⁰ + 5125309/2²⁴ ω + -1318839/2²⁵ ω² + -4192749/2²⁴ ω³  25541547/2²⁵ ω⁰ + 4095549/2²⁵ ω + 4061155/2²⁴ ω² + -6324583/2²⁵ ω³
```

Measure the execution time per gate in nano seconds.
```julia
julia> benchmark_compose(gates_r_pi_over_8)
86.59689922480621
```

Compute the error in operator norm

```julia
julia> using LinearAlgebra: opnorm;
julia> error = opnorm(RZ(big(pi)/8) - big(m)); @printf("%.3e", error)
5.017e-11
```
"""
function benchmark_compose end

end # module Benchmarking
