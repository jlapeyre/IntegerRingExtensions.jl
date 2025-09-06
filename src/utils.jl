module Utils

import LinearAlgebra
import ILog2: ilog2

export subscript, superscript

subscript(i::Integer) = _script(i, _sub_digit)
superscript(i::Integer) = _script(i, _super_digit)
superscript(x) = string("^", x)

function _script(i::Integer, _script_func)
    if i < 0
        error("Negative numbers are not supported for direct subscript conversion.")
    end
    # Get the digits of the integer in reverse order (e.g., 123 -> [3, 2, 1])
    digits_array = reverse(digits(i))
    # Convert each digit to its Unicode subscript character and join them
    return join(_script_func(d) for d in digits_array)
end

function _sub_digit(d::Integer)
    zerochar = 0x02080
    Char(zerochar + d)
end

function _super_digit(d::Integer)
    if d == 0
        Char(0x2070) # Superscript 0
    elseif d == 1
        Char(0x00B9) # Superscript 1
    elseif d == 2
        Char(0x00B2) # Superscript 2
    elseif d == 3
        Char(0x00B3) # Superscript 3
    elseif d in (4,5,6,7,8,9)
        Char(0x2070 + d) # Superscript for digits 4-9
    else
        throw(ArgumentError(lazy"Can't make superscript from digit '$d'"))
    end
end

const PRETTY = MIME"text/plain"

function _show_with_fieldnames(io::IO, obj)
    T = typeof(obj)
    print(io, T, "(")
    i = 0
    for field in fieldnames(T)
        i += 1
        i > 1 && print(io, ", ")
        print(io, field, " = ", getfield(obj, field))
    end
    print(io, ")")
end


function cpad(s::AbstractString, width::Integer, pad::AbstractString=" ")
    len = length(s)
    padlen = width - len
    if padlen <= 0
        return s
    end
    lpadlen = padlen ÷ 2 + padlen % 2
    rpadlen = padlen ÷ 2
    return lpad(s, len + lpadlen, pad) |> x -> rpad(x, width, pad)
end

function iszero_strong(x)
    iszero(x) === true
end

function isone_strong(x)
    isone(x) === true
end

function greater_than_strong(x::Number, y::Number)
    x > y
end

greater_than_strong(x, y) = false

"""
    lobit(z)

Return the position of the lowest set bit in `z`.

The rightmost bit (lsb) has position zero.
"""
function lobit(z)
    ilog2(z & -z)
end

"""
    small(x::BigFloat)::Float64

Convert `x` to `Float64` without getting `InexactError`.
"""
function small(x::BigFloat)
    Float64(trunc(x; sigdigits=17))
end

"""
    small(z::Complex{BigFloat})::Complex{Float64}

Convert `z` to `ComplexF64` without getting `InexactError`.
"""
function small(x::Complex{BigFloat})
    ComplexF64(trunc(x; sigdigits=17))
end

# This also exists in StatsBase

function countmap(itr)
    d = Dict{typeof(first(itr)), Int}()
    for x in itr
        d[x] = get(d, x, 0) + 1
    end
    d
end

# A copy of Base.power_by_squaring. Modified
# to not fail with Matrix2x2
function _power_by_squaring(x_, p::Integer; mul=*)
    x = x_
    if p == 1
        return x
    elseif p == 0
        return one(x)
    elseif p == 2
        return mul(x, x)
    elseif p < 0
        isone(x) && return x
        isone(-x) && return iseven(p) ? one(x) : x
        throw(ArgumentError(lazy"Can't raise object $x to power $p"))
    end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x = mul(x, x)
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x = mul(x, x)
        end
        y = mul(y, x)
    end
    return y
end

"""
    random_unitary(n::Integer)::Matrix{ComplexF64}

Best for `n>3`.
"""
function random_unitary(n::Integer)
    z = randn(ComplexF64, n, n)
    qr_fac = LinearAlgebra.qr(z)
    ph = LinearAlgebra.diagm([x / abs(x) for x in LinearAlgebra.diag(qr_fac.R)])
    return LinearAlgebra.mul!(z, qr_fac.Q, ph) # reuse storage in z for result
end

# This seems to work, but it's slow
"""
    random_special_unitary(n::Integer)

Return an `n`x`n` random special unitary matrix.

This is computed by scaling the eigenvalues of a Haar-random unitary matrix.
"""
function random_special_unitary(n::Integer)
    u = random_unitary(n)
    prod_eigs = LinearAlgebra.det(u)
    phase = prod_eigs^(-1/n)
    phase * u
end

end # module Utils
