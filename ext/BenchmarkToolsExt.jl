module BenchmarkToolsExt

using BenchmarkTools: @belapsed

import IntegerExtensions: compose, benchmark_compose

function benchmark_compose(gate_str)
    t = @belapsed compose($gate_str)
    return 10^9 * t / length(gate_str)
end

end # module BenchmarkToolsExt
