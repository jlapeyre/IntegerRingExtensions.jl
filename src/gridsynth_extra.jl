module GridSynthExtra

using ..CyclotomicRings: Domega
using ..GridSynth: GridSynthMatrix, GridSynthResults, parse_gridsynth_matrix
using ..Common: mul_root_two
using ..Matrices: Matrix2x2

function gridsynth_matrix_to_cyclic(matrix::GridSynthMatrix; power=true)
    domegas = map(x -> Domega(x...), matrix.data)
    if power
        domegas = map(x -> mul_root_two(x, -matrix.power), domegas)
    end
    Matrix2x2(domegas...,)
end

function gridsynth_matrix_to_cyclic(grid_results::GridSynthResults; power=true)
    gridsynth_matrix_to_cyclic(parse_gridsynth_matrix(grid_results); power=power)
end

end # module GridSynthExtra
