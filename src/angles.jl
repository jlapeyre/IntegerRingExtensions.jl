module Angles

"""
    Dar{T}

An angular unit that takes values on [-1, 1].
"""
struct Dar{T} <: Real
    function Dar(y)
        x = minus_one_to_one(y)
        new{typeof(x)}(x)
    end
    x::T
end

Base.:+(d::Dar, x::Real) = Dar(d.x + x)
Base.:-(d::Dar, x::Real) = Dar(d.x + x)
Base.:+(x::Real, d::Dar) = d + x
Base.:-(x::Real, d::Dar) = d - x
Base.:+(d1::Dar, d2::Dar) = Dar(d1.x + d2.x)
Base.:-(d1::Dar, d2::Dar) = Dar(d1.x - d2.x)

function zero_to_two(x)
    if x >= 0  # There must be a Base function call for this.
        (q, r) = divrem(x, 2)
        res =
            if iszero(q)
                r
            elseif iszero(r)
                2 * one(r)
            else
                r
            end
        (0 <= res <= 2) || error(lazy"sad outcome $ret")
        return res
    else
        n = round(x)
        y = x - 2 * n + 2
        y > 0 || error(lazy"I dont know what to do with $x")
        zero_to_two(y)
    end
end

function _minus_one_to_one(x)
    y = zero_to_two(x)
    if 0 <= y <= 1
        return y
    end
    return y - 2
end

function minus_one_to_one(x)
    y = _minus_one_to_one(x)
    (-1 <= y <= 1) || error(lazy"Sad case of $y")
    return y
end

# Probably want a Base.convert
function radtodar(theta)
    Dar(theta / pi)
end

function dartorad(dar::Dar)
    dar.x * pi
end

import Base: cos, sin, cis, tan

for func in (:cos, :sin, :cis, :tan)
    funcpi = Symbol(func, :pi)
    @eval $func(dar::Dar) = $funcpi(dar.x)
end

function Base.isapprox(a::Dar, b::Dar)
    isapprox(a.x, b.x)
end

function sadcheck(x, y)
    if  !isapprox(Dar(x + y),  (Dar(x) + Dar(y)))
#        @show (x, y)
        return true
    end
    return false
end

mrand() = 20 * rand()
#mrand() = 10 * rand()

sadcheck() = sadcheck(mrand(), mrand())

function runsad(N)
    sum(sadcheck() for _ in 1:N)
end

end # module Angles

