# Test that shift_minus_one_to_one(x)
# gives the same results as the elaborate function in angles.jl

function shift_minus_one_to_one(x)
    (f, w) = modf(x) # fractional, whole
    iszero(f) && return iseven(x) ? zero(x) : one(x)
    isodd(w) && (w += sign(x))
    return x - w
end

function test_shift_funcs(x)
    !isapprox(shift_minus_one_to_one(x), minus_one_to_one(x))
end

mrand() = 20*(rand() - 1/2)

function countbad(N)
    sum(_ -> test_shift_funcs(mrand()), 1:N)
end

function printbad(N)
    for _ in 1:N
        x = mrand()
        if test_shift_funcs(x)
            @show x
        end
    end
end
