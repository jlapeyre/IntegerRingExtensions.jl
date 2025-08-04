module GridSynth

using ..Utils: PRETTY

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

# function GridSynthOpts(theta::String;
#                        epsilon::Float64=1e-10, uptophase::Bool=false,
#                        effort::Int=25, hexoutput=false,
#                        stats::Bool=true, latex::Bool=false,
#                        seed::Int=0)
#     GridSynthOpts(theta, epsilon, uptophase, effort, hexoutput, stats, latex, seed)
# end

function makecommand(opts::GridSynthOpts)
    (;theta, epsilon, uptophase, effort, hexoutput, stats, latex, seed) = opts
    command = `gridsynth $theta -e $epsilon -f $effort`
    uptophase && push!(command.exec, "-p")
    hexoutput && push!(command.exec, "-h")
    stats && push!(command.exec, "-s")
    latex && push!(command.exec, "-l")
    seed > 0 && push!(command.exec, "-s", seed)
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
    # _s(obj) = show(io, PRETTY(), obj)
    # props = propertynames(gr)
    # println(io, typeof(gr), "(")
    # for (i, name) in enumerate(props)
    #     print(io, "  ", name, "=")
    #     prop = getproperty(gr, name)
    #     _s(prop)
    #     print(io, "\n")
    # end
    # print(io, ")")
end

_show(io::IO, obj, n) = show(io, PRETTY(), obj)

function _show(io::IO, gr::Union{GridSynthResults, GridSynthOpts}, n::Int)
#    _s(obj) = _show(io, PRETTY(), obj)
    _s(obj) = _show(io, obj, n+2)
    props = propertynames(gr)
    println(io, typeof(gr), "(")
    for (i, name) in enumerate(props)
        for i in 1:n
            print(io, " ")
        end
        print(io, name, "=")
        prop = getproperty(gr, name)
        _s(prop)
        print(io, "\n")
    end
    for i in 1:(n-1)
        print(io, " ")
    end
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
    theta = String(lspl(_theta))
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

function gridsynth(;kwargs...)
    opts = GridSynthOpts(;kwargs...)
    run_gridsynth(opts)
end

end # module GridSynth
