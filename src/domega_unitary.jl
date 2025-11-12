module DOmegaUnitaries

import ..RootOnes: Omega, omega
import ..CyclotomicRings: DOmega
import ..Singletons: InvRootTwo

export  DOmegaUnitary, mul_by_T_from_left, mul_by_T_inv_from_left, mul_by_H_and_T_power_from_left,
    mul_by_T_power_from_left, mul_by_S_from_left, mul_by_S_power_from_left,
    mul_by_W

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

# Multiply the bottom row of m by fac
function _mul_bottom_row(m::DOmegaUnitary, fac::Omega)
    DOmegaUnitary(m.u, fac * m.t, fac * m.omega_pow)
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
mul_by_T_from_left(m::DOmegaUnitary) = _mul_bottom_row(m, omega)
mul_by_T_inv_from_left(m::DOmegaUnitary) = _mul_bottom_row(m, inv(omega))
mul_by_T_power_from_left(m::DOmegaUnitary, n::Integer) = _mul_bottom_row(m, omega^n)
mul_by_S_from_left(m::DOmegaUnitary) = _mul_bottom_row(m, omega^2)
mul_by_S_power_from_left(m::DOmegaUnitary, n::Integer) = _mul_bottom_row(m, omega^(2*n))

# TODO: reduce elements, probably
function mul_by_H_from_left(m::DOmegaUnitary)
    (u, t, omega_pow) = (m.u, m.t, m.omega_pow)
    DOmegaUnitary(InvRootTwo * (u + t), InvRootTwo * (u - t), -omega_pow)
end

function mul_by_H_and_T_power_from_left(m::DOmegaUnitary, n::Integer)
    mul_by_H_from_left(mul_by_T_power_from_left(m, n))
end

function mul_by_X_from_left(m::DOmegaUnitary)
    DOmegaUnitary(m.t, m.u, -m.omega_pow)
end

# Multiply by global factor of ω
function mul_by_W(m::DOmegaUnitary)
    DOmegaUnitary(omega * m.u, omega * m.t, omega^2 * m.omega_pow)
end

# Multiply by global factor of ωⁿ
function mul_by_W_power(m::DOmegaUnitary)
    DOmegaUnitary(omega * m.u, omega * m.t, omega^2 * m.omega_pow)
end

end # module DOmegaUnitaries
