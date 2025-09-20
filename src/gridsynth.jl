module GridSynth

using DispatchDoctor: @stable

using ..Utils: PRETTY, superscript
import ..Utils: countmap
import ..Compose: Compose, compose, compose_scale, rotation_error

# Parameters disables the default constructor :(
# Also does not allow forwarding keywords args with (;kwargs...)
# Base.@kwdef supports both of these.
#
# using Parameters: @with_kw, Parameters
# Parameters.@with_kw struct GridSynthOpts

@kwdef struct GridSynthOpts
    theta::String
    epsilon::Float64 = 1e-10 # precision
    uptophase::Bool = false
    effort::Int = 25
    hexoutput::Bool = false
    stats::Bool = true
    latex::Bool = false
    seed::Int = 0
end

function makecommand(opts::GridSynthOpts)
    (;theta, epsilon, uptophase, effort, hexoutput, stats, latex, seed) = opts
    command = `gridsynth -e $epsilon -f $effort`
    uptophase && push!(command.exec, "-p")
    hexoutput && push!(command.exec, "-h")
    stats && push!(command.exec, "-s")
    latex && push!(command.exec, "-l")
    seed > 0 && push!(command.exec, "-s", seed)
    push!(command.exec, "--")
    push!(command.exec, "$theta")
    return command
end

@kwdef struct GridSynthResults
    gates::String
    seed::Tuple{Int, Int}
    tcount::Int
    lb_tcount::Int
    theta::String
    epsilon::String
    matrix::String
    error::String
    runtime::Float64 # in seconds
    candidates::@NamedTuple{failed::Int, timedout::Int, succeeded::Int}
    timepercandidate::Float64
    opts::GridSynthOpts
end

function Base.show(io::IO, ::PRETTY, gr::Union{GridSynthResults, GridSynthOpts})
    _show(io, gr, 1)
end

_show(io::IO, obj, n) = show(io, PRETTY(), obj)

function _show(io::IO, gr::Union{GridSynthResults, GridSynthOpts}, n::Int)
    _s(obj) = _show(io, obj, n+2)
    props = propertynames(gr)
    println(io, typeof(gr), "(")
    for (i, name) in enumerate(props)
        print(io, " "^n)
        print(io, name, "=")
        prop = getproperty(gr, name)
        _s(prop)
        print(io, "\n")
    end
    print(io, " "^(n-1))
    print(io, ")")
end


run_gridsynth(command::Cmd) = split(read(command, String), "\n")

function run_gridsynth(opts::GridSynthOpts)
    lines = run_gridsynth(makecommand(opts))
    store_results(lines, opts)
end

function store_results(lines, opts)
    (gates, _seed, _tcount, _lb, _theta, _epsilon, _matrix,
     _error, _runtime, _candidates, _timecandidate) = lines
    # GridSynthResults(
    #     lines
    # )
    spl = s -> split(s, r"\s+")
    ps = s -> parse(Int, s)
    lspl = s -> last(spl(s))

    seed = Tuple((ps(x) for x in spl(_seed)[3:4]))
    tcount = ps(lspl(_tcount))
    lowerbound = ps(lspl(_lb))
    theta = String(replace(_theta, r"Theta:\s+" => ""))
    epsilon = String(lspl(_epsilon))

    _matrix_res  = collect(split(_matrix, r":"))
    length(_matrix_res) == 2 || throw(ErrorException(lazy"error parsing matrix"))
    matrix = String(lstrip(last(_matrix_res)))
    error = String(lspl(_error))

    _c1 = match(r"\((.+)\)", _candidates).captures[1]
    cands = ps.(String.(first.(split.(split(_c1, r",\s"), r"\s"))))
    candidates = (failed=cands[1], timedout=cands[2], succeeded=cands[3])

    timecandidate = last(spl(_timecandidate))
    timecandidates = parse(Float64, timecandidate[1:end-1])

    runtime = last(spl(_runtime))
    runtime = parse(Float64, runtime[1:end-1])

    gates = String(gates)

    lb_tcount=lowerbound
    timepercandidate = timecandidates
    rtup = (gates, seed, tcount, lowerbound, theta, epsilon, matrix, error, runtime, candidates, timecandidates)

    # Note we are using magic keyword assignment. The names of the arguments match the fields in the structs.
    # The arguments here can be supplied in any order.
    GridSynthResults(;
                     gates,
                     seed,
                     tcount,
                     lb_tcount,
                     theta,
                     epsilon,
                     matrix,
                     error,
                     runtime,
                     candidates,
                     timepercandidate,
                     opts
    )
end

# TODO: Convert float inputs to strings
function gridsynth(;kwargs...)
    # @show typeof(kwargs)
    # @show kwargs
    # if haskey(kwargs, :theta)
    #     kwargs[:theta] = string(kwargs[:theta])
    # end
#    kwargs = Tuple(
    opts = GridSynthOpts(;kwargs...)
    run_gridsynth(opts)
end

###
### Lossless or arbitrary precision conversion of numbers as strings to number types.
###

function stringtonum(gr::GridSynthResults)
    return (theta = biggennum(gr.theta), epsilon = biggennum(gr.epsilon), error = biggennum(gr.error))
#    return (theta = stringtonum(gr.theta), epsilon = stringtonum(gr.epsilon), error = stringtonum(gr.error))
end

using MacroTools: postwalk, @capture

# Change big(b)^big(n) to big(b)^n  . `n` can be an expr like `-n`
# If n is a negative literal, evaluation will fail if it's made big.
function unbigpow(expr::Expr)
    postwalk(x -> @capture(x, ^(big(b_), big(n_))) ? :(^(big($b), $n)) : x, expr)
end

biggennum(str::AbstractString) = eval(biggenex(str))

function biggenex(str::AbstractString)
    # First modify expressions output by gridsynth to be valid Julia expressions
    str = replace(str, r"(\*\*)" => "^")
    str = replace(str, r"sqrt\s+([^\s]+)" => s"sqrt(\1)") # fix square root

    # Next use MacroTools to rewrite the expression.
    biggenex(Meta.parse(str))
end

function biggenex(expr::Expr)
    res = postwalk(expr) do x
        if isa(x, Integer)
            return :(big($x))
        end
        if x === :pi
            return :(big(pi))
        end
        if isa(x, Float64)
            xs = string(x)
            return :(@big_str($xs))
        end
        return x
    end
    unbigpow(res)
end

function compose(gr::GridSynthResults; chunklen=nothing)
    return isnothing(chunklen) ? compose(gr.gates) : compose(gr.gates; chunklen)
end

function compose_scale(gr::GridSynthResults; chunklen=nothing)
    return isnothing(chunklen) ? compose_scale(gr.gates) : compose_scale(gr.gates; chunklen)
end

function rotation_error(gr::GridSynthResults)
    rotation_error(compose(gr), biggennum(gr.theta))
end

function countmap(gr::GridSynthResults)
    countmap(gr.gates)
end

struct GridSynthMatrix{T}
    data::Matrix{T}
    power::Int
end

function Base.show(io::IO, ::PRETTY, m::GridSynthMatrix)
    print(io, typeof(m), ": 1/√2", superscript(m.power),  " x\n")
    show(io, PRETTY(), m.data)
end

"""
    gridsynth_matrix(str::AbstractString, ::Type{T}=BigInt) where {T}

Parse the string representation of a 2x2 matrix over `𝔻[ω] = ℤ[1/√2, i]` and return
a numeric representation.

The representation is a `Matrix{Vector{BigInt}}` and a power of one over root two.
The power is of type `Int64`. Each element is a `Vector` of coefficients of
an element of 𝔻[ω]

The returned value is a `GridSynthMatrix{T}`.
"""
function gridsynth_matrix(str::AbstractString, ::Type{CoeffT}=BigInt) where {CoeffT}
    rhr = r"roothalf\^(\d+)\s\*\s"
    rootmatch = match(rhr, str)
    if !isnothing(rootmatch)
        roothalf = parse(Int, only(rootmatch.captures))
    else
        roothalf = 0
    end
    str = replace(str, rhr => "")
    str = replace(str, r"matrix\s+" => "")
    rows = collect(eachmatch(r"(\[Omega[^\]]+\])", str))
    length(rows) == 2 || error(lazy"Expected two rows")
    row1 = only(rows[1].captures)
    row2 = only(rows[2].captures)
    # Reverse because S+R choose opposite order for coefficients
    getcoeffs(el) = reverse!([parse(CoeffT, String(only(x.captures))) for x in eachmatch(r"(-?\d+)", el)])
    matels = [split(row1, ",")..., split(row2, ",")...]
    matrix = permutedims(collect(reshape(map(getcoeffs, matels), (2,2))))
    return GridSynthMatrix(matrix, roothalf)
end

function gridsynth_matrix(grid_results::GridSynthResults, ::Type{CoeffT}=BigInt) where {CoeffT}
    gridsynth_matrix(grid_results.matrix, CoeffT)
end


end # module GridSynth
