@stable module GridSynthExtra

using ..CyclotomicRings: DOmega
using ..GridSynth: GridSynthMatrix, GridSynthResults, gridsynth_matrix
using ..Common: mul_root_two
using ..Matrices2x2: Matrix2x2

function gridsynth_matrix_to_cyclic(matrix::GridSynthMatrix; power=true)
    domegas = map(x -> DOmega(x...), matrix.data)
    if power
        domegas = map(x -> mul_root_two(x, -matrix.power), domegas)
    end
    Matrix2x2(domegas...,)
end

function gridsynth_matrix_to_cyclic(grid_results::GridSynthResults; power=true)
    gridsynth_matrix_to_cyclic(gridsynth_matrix(grid_results); power=power)
end

end # module GridSynthExtra
