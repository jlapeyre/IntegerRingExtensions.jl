module DOmegaUnitaries

import ..RootOnes: Omega, omega
import ..CyclotomicRings: DOmega

export mul_by_T_from_left, mul_by_T_inv_from_left, mul_by_H_and_T_power_from_left

const pretty = MIME"text/plain"

# Represent a 2x2 unitary matrix over `𝔻[ω] = ℤ[1/√2, i]`:
# u  -t^* ωᵏ
# t   u^* ωᵏ
struct DOmegaUnitary{T}
    u::DOmega{T}
    t::DOmega{T}
    omega_pow::Omega # ωᵏ
end

function Base.:(==)(m1::DOmegaUnitary, m2::DOmegaUnitary)
    m1.u == m2.u && m1.t == m2.t && m1.omega_pow == m2.omega_pow
end

function Base.show(io::IO, ::pretty,  m::DOmegaUnitary)
    (u, t) = (m.u, m.t)
    spcs = " "^5
    show(io, pretty(), u)
    print(io, spcs)
    show(io, pretty(), -adjoint(t))
    println(io)
    show(io, pretty(), t)
    print(io, spcs)
    show(io, pretty(), adjoint(u))
    println(io)
    nothing
end

# Mulitplying by T from the left multiplies the bottom row by omega.
# Before multiplying we have:
# u  -t^* ωᵏ
# t   u^* ωᵏ
# Instantiating with u, ωt, ωᵏ⁺¹ we have
# u    -t^* ω^* ωᵏ⁺¹     u    -t^* ωᵏ
#                    ==
# tω   u^* ωᵏ⁺¹          tω   u^* ωᵏ⁺¹
# This has the effect of multiplying the bottom row by ω, as desired.
function mul_by_T_from_left(m::DOmegaUnitary)
    DOmegaUnitary(m.u, omega * m.t, omega * m.omega_pow)
end

function mul_by_T_inv_from_left(m::DOmegaUnitary)
    DOmegaUnitary(m.u, inv(omega) * m.t, inv(omega) * m.omega_pow)
end

function mul_by_T_power_from_left(m::DOmegaUnitary, n::Integer)
    DOmegaUnitary(m.u, omega^n * m.t, omega^n * m.omega_pow)
end

end # module DOmegaUnitaries
