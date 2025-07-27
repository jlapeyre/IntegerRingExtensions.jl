module FunctionSecants

export secant_slope, secant

#const _SPECIALIZED_TYPES_REAL = Union{typeof(log)}
#const _SPECIALIZED_TYPES = Union{typeof(cos), typeof(sin), typeof(exp), typeof(cis), typeof(cispi)}
const _SPECIALIZED_TYPES = Union{typeof(cos), typeof(sin), typeof(exp), typeof(cis), typeof(cispi), typeof(log), typeof(sqrt)}

_d_small(d) = abs(d) < sqrt(eps(real(typeof(d))))

function secant_slope(func::typeof(sqrt), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    if _d_small(d)
        inv(sqrt(w)) * evalpoly(d/w, (1, 0, 1//8, 0, 7//128, 0, 33//1024, 0, 715//32768, 0, 4199//262144))
    else
        (func(x) - func(y))/ (x-y)
    end
end

function _atanh_div_x(x)
    abs(x) < eps(real(typeof(x))) && return one(x)
    atanh(x)/x
end

# Some inputs give good results, others very wrong results.
function secant_slope(::typeof(log), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return inv(w) * _atanh_div_x(d / w)
end

function _get_wd(x::T, y::T) where {T}
    ((x + y) / 2, (x - y) / 2)
end

function _sinh_div_x(x)
    abs2(x) < eps(real(typeof(x)))^2 && return one(x)
    sinh(x)/x
end

function secant_slope(::typeof(exp), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return exp(w) * _sinh_div_x(d)
end

function secant_slope(::typeof(cispi), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return im * pi * cispi(w) * sinc(d)
end

function secant_slope(::typeof(cis), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return im * cis(w) * sinc(d / pi)
end

function secant_slope(::typeof(sin), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return cos(w) * sinc(d / pi)
end

function secant_slope(::typeof(cos), x::T, y::T) where {T <: Complex}
    (w, d) = _get_wd(x, y)
    return -sin(w) * sinc(d / pi)
end

function secant_slope(::typeof(cos), x::T, y::T) where {T <: Real}
    (w, d) = _get_wd(x, y)
    return -sin(w) * sinc(d / pi)
end

function secant(f::_SPECIALIZED_TYPES, x::T, y::T) where {T <: Number}
    return (f(x), f(y), secant_slope(f, x, y))
end

# function secant(f::_SPECIALIZED_TYPES_REAL, x::T, y::T) where {T <: Real}
#     return (f(x), f(y), secant_slope(f, x, y))
# end

function _central_diff(f, x, z, h=1e-6)
    w = (x + z)/2
    (f(w+h) - f(w-h)) / (2*h)
end

fakeabs(x::Number) = abs(x)
fakeabs(z::Complex) = abs(real(z)) + abs(imag(z))

"""
    secant(f::F, x::T, y::T) where {F, T <: Number}

Return `(f(x), f(y), slope) where `slope = (f(x) - f(y)) / (x - y)`.

`slope` is determined by the function `secant_slope`.
"""
function secant(f::F, x::T, y::T) where {F, T <: Number}
    fx = f(x)
    fy = f(y)
    dl = x - y
    if fakeabs(dl) < sqrt(eps(real(typeof(x))))
        sl = _central_diff(f, x, y)
    else
        sl = (fx - fy) / (x - y)
    end
    (fx, fy, sl)
end

"""
    secant_slope(f::F, x::T, y::T) where {F, T <: Number}

Return `(f(x) - f(y)) / (x - y)`.

Arguments such that `abs(x - y)` is small, or vanishes, are supported to varying degrees of accuracy.
For several functions `f`, a special implementation gives high accuracy in this case.

More generic, high accuracy algorithms exist, but are not implemented here.
"""
function secant_slope(f::F, x::T, y::T) where {F, T <: Number}
    (fx, fy, sl) = secant(f, x, y)
    return sl
end

end  # module FunctionSecants
