module Gates

# Better to use relative path. If I could find out how
using ..IntegerExtensions.Matrices: Matrix2x2
using ..IntegerExtensions: imaginary, sqrt_imaginary, one_over_root_two
using ..IntegerExtensions.Rings: Domega, canonical
using Nemo: ZZ, ZZRingElem

function Zgate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    Matrix2x2(o, z, z, -o)
end

function Sgate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    img = imaginary(T)
    Matrix2x2(o, z, z, img)
end

function Tgate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    sqrt_img = sqrt_imaginary(T)
    Matrix2x2(o, z, z, sqrt_img)
end

function Hgate(::Type{T}) where T
    inv_half = one_over_root_two(T)
    Matrix2x2(inv_half, inv_half, inv_half, -inv_half)
end

function Xgate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    Matrix2x2(z, o, o, z)
end

function Igate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    Matrix2x2(o, z, z, o)
end

function gate_map(::Type{T}) where T
   return Dict(
        :H => Hgate(T),
        :S => Sgate(T),
        :Z => Zgate(T),
        :X => Xgate(T),
        :T => Tgate(T),
        :I => Igate(T)
    )
end

const GATE_MAP_BIG = gate_map(Domega{BigInt})
const GATE_MAP_INT = gate_map(Domega{Int})
const GATE_MAP_ZZ = gate_map(Domega{ZZRingElem})

function compose(gates::AbstractString, gmap=GATE_MAP_BIG)
    result = gmap[:I]
    for gate in Iterators.reverse(gates)
        result = gmap[Symbol(gate)] * result
    end
    return result
end

end # module Gates
