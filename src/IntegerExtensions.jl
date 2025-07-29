module IntegerExtensions
import Reexport

include("rings.jl")
Reexport.@reexport import .Rings: QuadraticRing, isunit, RootOne, DyadicFraction, CyclotomicRing, RootOneA,
    RootOne8, QuadraticRing2, canonical

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, power

end # module IntegerExtensions
