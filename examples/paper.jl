using IntegerExtensions

const delta = Domega(1, 1, 0 , 0)

const idelta = mul_one_over_root_two(Domega(0, 1, -1, 0))
@assert isone(delta * idelta)

const cdelta = complex(delta)
const cidelta = complex(idelta)

# Can represent with Zomega. But we need to manipulate delta as a Domega
const delta_z = Zomega(1, 1, 0 , 0)


const lambda = Zroot2(1, 1)

#const delta = CyclotomicRing{4, Int}(0, 0, 1, 1)
#const lambda = QuadraticRing2(1, 1)

nothing;
