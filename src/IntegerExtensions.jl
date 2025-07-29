module IntegerExtensions
import Reexport

export Domega

include("rings.jl")
Reexport.@reexport import .Rings: QuadraticRing, isunit, RootOne, DyadicFraction, CyclotomicRing, RootOneA,
    RootOne8, QuadraticRing2, canonical, Domega

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, power

include("gates.jl")
Reexport.@reexport import .Gates: Zgate

end # module IntegerExtensions
