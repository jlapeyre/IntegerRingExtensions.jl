module GateMatrix

using ..Gates: Gate1
using ..Compose: compose
using ..Matrices2x2: Matrix2x2
using ..RootOnes: Omega

for g in (:X, :Y, :Z, :S, :T)
    gn = Symbol(g, :F64)
    @eval const $gn = Matrix2x2{ComplexF64}(Gate1{$(QuoteNode(g))})
end

# Not really isclifford
# function isclifford(m)
#     for g in (XF64, YF64, ZF64)
#         show(stdout, PRETTY(), zchop((m' * g * m') * (1+im)))
#         println("\n")
#     end
# end

const CLIFFORD_STR = (
    :I,
    :X,
    :Y,
    :Z,
    :S,
    :SX,
    :SY,
    :SZ,
    :H,
    :HX,
    :HY,
    :HZ,
    :HS,
    :HSX,
    :HSY,
    :HSZ,
    :SH,
    :SHX,
    :SHY,
    :SHZ,
    :SHS,
    :SHSX,
    :SHSY,
    :SHSZ)

const CLIFFORD_DOMEGA = map(s -> compose(string(s)), CLIFFORD_STR)

const CLIFFORD_DOMEGA_MAP = Dict{Symbol, typeof(CLIFFORD_DOMEGA[1])}()

for (s, m) in zip(CLIFFORD_STR, CLIFFORD_DOMEGA)
    CLIFFORD_DOMEGA_MAP[s] = m
end

function find_clifford(m)
    for c in CLIFFORD_DOMEGA
        for k in 0:7
            if Omega(k) * c == m
                return (c, k)
            end
        end
    end
    return nothing
end

end # module GateMatrix
