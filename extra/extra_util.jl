# I don't want to depend on ZChop. Could make this an extension of something.
# For now, just  load this with `include`

using ZChop

using IntegerExtensions: Matrix2x2

function ZChop.zchop(m::Matrix2x2)
    (a, b, c, d) = m.data
    Matrix2x2(zchop(a), zchop(b), zchop(c), zchop(d))
end
