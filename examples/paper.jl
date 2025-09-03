using IntegerExtensions

const delta = DOmega(1, 1, 0 , 0)

const idelta = mul_one_over_root_two(DOmega(0, 1, -1, 0))
@assert isone(delta * idelta)

const cdelta = complex(delta)
const cidelta = complex(idelta)

# Can represent with ZOmega. But we need to manipulate delta as a DOmega
const delta_z = ZOmega(1, 1, 0 , 0)

const lambda = Zroot2(1, 1)

@assert abs(float(delta^2)) == 2 + sqrt(2)

nothing;
