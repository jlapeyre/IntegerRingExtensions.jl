module Utils

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


end # module Utils
