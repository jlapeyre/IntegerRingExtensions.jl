module Gates

# Better to use relative path. If I could find out how
using ..Matrices: Matrix2x2
using ..Common: imaginary, sqrt_imaginary, one_over_root_two, canonical
using ..CyclotomicRings: Domega
using ..QuadraticRings: Droot2

# using Nemo: ZZ, ZZRingElem

function Igate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    Matrix2x2(o, z, z, o)
end

function Xgate(::Type{T}) where T
    o = one(T)
    z = zero(T)
    Matrix2x2(z, o, o, z)
end

function Ygate(::Type{T}) where T
    z = zero(T)
    img = imaginary(T)
    Matrix2x2(z, img, -img, z)
end

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


function Wgate(::Type{T}) where T
    z = zero(T)
    sqrt_img = sqrt_imaginary(T)
    Matrix2x2(sqrt_img, z, z, sqrt_img)
end

function gate_map(::Type{T}) where T
   return Dict(
       :H => Hgate(T),
       :S => Sgate(T),
       :Z => Zgate(T),
       :X => Xgate(T),
       :T => Tgate(T),
       :I => Igate(T),
       :W => Wgate(T)
    )
end

const GATE_MAP_INT = gate_map(Domega{Int})
const GATE_MAP_BIG_INT = gate_map(Domega{BigInt})
const GATE_MAP_INT128 = gate_map(Domega{Int128})
const GATE_MAP_BIG_FLOAT = gate_map(BigFloat)

const GATE_MAP_GOOD = Dict(
    :H => Hgate(Droot2{Int, Int}),
    :S => Sgate(Int),
    :X => Xgate(Int),
    :T => Tgate(Domega{Int}),
    :I => Igate(Domega{Int}),
    :W => Wgate(Domega{Int}),
)

# Careful this may segfault when using precompiled because it
# calls a dynamically linked C library.
# const GATE_MAP_ZZ = gate_map(Domega{ZZRingElem})

"""
    compose(gates::AbstractString, gmap=GATE_MAP_BIG_INT; reduce_fractions=true)

Compute composition of the gates in `gates`.

`gmap` is a map from `Symbol`s to matrices. If `reduce_fractions` is `true`
then reduce fractions in `DyadicFraction`s after each matrix multiplication.

`reduce_fractions` reduces the maximum values of intermediate numbers allowing computation
of longer compositions with smaller data types.
"""
function compose(gates_in::AbstractString)
#function compose(gates_in::AbstractString; reduce_fractions=true)
    chunklen = 300
#    reduce_func = reduce_fractions ? canonical : identity
    chunks = reverse(chunkstring(reverse(gates_in), chunklen))
    mats = [map(Domega{BigInt}, compose_one(chunk, GATE_MAP_INT)) for chunk in chunks]
    canonical(prod(mats))
end

"""
    chunkstring(s, chunklen=5)

Return an array of chunks of string `s`, each of length `chunklen`,
except the last one, which may be shorter. Satisfies
```
join(chunkstring(s)) == s
```
"""
function chunkstring(s, chunklen=5)
    strs = String[]
    i = 1
    while true
        cend = i + chunklen - 1
        cend = cend < length(s) ? cend : length(s)
        push!(strs, s[i:cend])
        cend == length(s) && break
        i = cend + 1
    end
    return strs
end

# Assume gates has already been reversed!
"""
    compose_one(gates::AbstractString, gmap=GATE_MAP_BIG_INT; reduce_fractions=true)

Compose gates in string `gates`. This is meant to compute the composition for a chunk
of a longer string of gates. The chunking and recombining is done by `compose`.
The reason we compute by chunks is that the chunks are small enough that computation
can be done with 64-bit (i.e. fast) integers. Then the resulting matrices are converted to
`BigInt`, and a final composition of these matrices is peformed.

Compute composition of the gates in `gates`. The composition will be from
*right to left*. In particular, you need to reverse the output string from
`gridsynth` before calling `compose`.

`gmap` is a map from `Symbol`s to matrices. If `reduce_fractions` is `true`
then reduce fractions in `DyadicFraction`s after each matrix multiplication.

`reduce_fractions` reduces the maximum values of intermediate numbers allowing computation
of longer compositions with smaller data types.
"""
function compose_one(gates::AbstractString, gmap=GATE_MAP_BIG_INT; reduce_fractions=true)
    result = gmap[:I]
    reduce_func = reduce_fractions ? canonical : identity
    for gate in gates
        result = reduce_func(gmap[Symbol(gate)] * result)
    end
    return result
end

"""
    RZ(theta)

This is the Z rotation matrix that is synthesized by `gridsynth` as
well as algorithms in other papers.
"""
function RZ(theta)
    z = zero(theta)
    t2 = theta/2
    Matrix2x2(cis(-t2), z, z, cis(t2))
end

end # module Gates
