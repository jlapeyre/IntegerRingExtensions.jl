module RingMatrices

import ..CyclotomicRings: Domega
import ..Matrices2x2: Matrix2x2
import ..RootOnes: RootOne

# The modules for rings, groups, and Matrices are independent.
# This module mixes them.

function Base.:*(r::RootOne{8}, m::Matrix2x2{<:Domega})
    map(x -> r * x, m)
end

end #module RingMatrices
