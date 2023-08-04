using MPI
MPI.Init()

"""Wrappers around the logging macros to ensure they only run for a single process.."""

const comm::MPI.Comm = MPI.COMM_WORLD

macro mpidebug(exs...)
    return :(MPI.Comm_rank(comm) > 0 ? nothing : @debug $(esc.(exs)...))
end

macro mpiinfo(exs...)
    return :(MPI.Comm_rank(comm) > 0 ? nothing : @info $(esc.(exs)...))
end

macro mpiwarn(exs...)
    return :(MPI.Comm_rank(comm) > 0 ? nothing : @warn $(esc.(exs)...))
end

macro mpierror(exs...)
    return :(MPI.Comm_rank(comm) > 0 ? nothing : @error $(esc.(exs)...))
end

# This was so troublesome until I found https://discourse.julialang.org/t/loggingextras-jl-how-to-add-custom-log-macro/64679/7
function make_kwargs(ex)
    if ex isa Expr && ex.head === :(=) ||
       ex isa Expr && ex.head === :...
        return ex
    else
        Expr(:(=), ex, esc(ex))
    end
end

function get_rank_msg()::String
    return "[Process $(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))]"
end

macro mpirankeddebug(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Debug
    message = get_rank_msg() * " " * string(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpirankedinfo(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Info
    message = get_rank_msg() * " " * string(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpirankedwarn(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Warn
    message = get_rank_msg() * " " * string(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpirankederror(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Error
    message = get_rank_msg() * " " * string(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

