using ..CyclotomicRings: DOmega, ZOmega, residue

const BIT_SHIFT = (0, 0, 1, 0, 2, 0, 1, 3, 3, 3, 0, 2, 2, 1, 0, 0)
const BIT_COUNT = (0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4)

const _T_POWER_AND_H = ((:H,), (:T, :H), (:S, :H), (:T, :S, :H))

function t_power_and_h(m)
    _T_POWER_AND_H(m+1)
end

function _reduce_denomexp(unitary::Unitary2x2)
    residue_u = residue(unitary_u(unitary))
    residue_t = residue(unitary_t(unitary))
    residue_squared_u = residue(residue_u * residue_u')
    m = BIT_SHIFT[residue_t] - BIT_SHIFT[residue_u]
    m < 0 && (m += 4)
    if iszero(residue_squared_u)
        unitary = unitary.mul_by_H_and_T_power_from_left(0).renew_denomexp(
            unitary.k - 1
        )
        return t_power_and_h(0), unitary
    end
    if residue_squared_z == 0b1010
        unitary = unitary.mul_by_H_and_T_power_from_left(-m).renew_denomexp(
            unitary.k - 1
        )
        return t_power_and_h(m), unitary
    end
end
