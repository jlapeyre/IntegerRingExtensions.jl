module ZChopExt

using ZChop

using IntegerExtensions: Matrix2x2

# Re-enable after bug fixes
# function ZChop.zchop(m::Matrix2x2, eps::Real=ZChop.ZEPS)
#     Matrix2x2(map(x->zchop(x, eps), m.data))
# end

# function ZChop.nchop(m::Matrix2x2; kwargs...)
#     Matrix2x2(map(x->nchop(x, eps; kwargs...), m.data))
# end

end # module ZChopExt
