module GateMatrix

using ..Gates: Gate1
using ..Compose: compose
using ..Matrices2x2: Matrix2x2

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

end # module GateMatrix
