module Singletons

import Base: show, inv
import ..Utils: PRETTY

struct Root2T end

"""
    Root2

Represents the square root of two.
"""
const Root2 = Root2T()

show(io::IO, ::PRETTY, ::Root2T) = print(io, "√2")

struct InvRoot2T end

"""
    InvRoot2

Represents the reciprocoal of the square root of two.
"""
const InvRoot2 = InvRoot2T()

show(io::IO, ::PRETTY, ::InvRoot2T) = print(io, "√2⁻¹")

inv(::InvRoot2T) = Root2
inv(::Root2T) = InvRoot2

struct RootIT end

"""
    RootI

Represents the principal square root of the imaginary unit.
This is also the principal eight root of one.
"""
const RootI = RootIT()

show(io::IO, ::PRETTY, ::RootIT) = print(io, "√𝕚")

end # module Singleons
