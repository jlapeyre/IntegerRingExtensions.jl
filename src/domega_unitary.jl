module DOmegaUnitaries

import ..CyclotomicRings: DOmega

struct DOmegaUnitary{T}
    n::Int
    u::DOmega{T}
    t::DOmega{T}
end

end # module DOmegaUnitaries
