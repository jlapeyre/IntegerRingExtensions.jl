module RingMatrices

import ..CyclotomicRings: DOmega
import ..Matrices2x2: Matrix2x2
import ..RootOnes: RootOne

# The modules for rings, groups, and Matrices are independent.
# This module mixes them.

function Base.:*(r::RootOne{8}, m::Matrix2x2{<:DOmega})
    map(x -> r * x, m)
end

Base.:*(m::Matrix2x2{<:DOmega}, r::RootOne{8}) = r * m

function compute_phase_factor(m1::Matrix2x2{<:DOmega}, m2::Matrix2x2{<:DOmega})
    for i in 0:7
        RootOne{8}(i) * m1 == m2 && return i
    end
    return -1
end

end #module RingMatrices
