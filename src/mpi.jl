using MPI
MPI.Init()

"""Wrappers around the logging macros to ensure they only run for a single process.."""

const comm::MPI.Comm = MPI.COMM_WORLD

# This was so troublesome until I found https://discourse.julialang.org/t/loggingextras-jl-how-to-add-custom-log-macro/64679/7
function make_kwargs(ex)
    if ex isa Expr && ex.head === :(=) ||
       ex isa Expr && ex.head === :...
        return ex
    else
        Expr(:(=), ex, esc(ex))
    end
end

macro mpidebug(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Debug
    message = esc(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    MPI.Comm_rank(comm) > 0 ? nothing : Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpiinfo(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Info
    message = esc(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    MPI.Comm_rank(comm) > 0 ? nothing : Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpiwarn(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Warn
    message = esc(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    MPI.Comm_rank(comm) > 0 ? nothing : Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
end

macro mpierror(exs...)
    location = LineNumberNode(__source__.line, __source__.file)
    level = Base.CoreLogging.Error
    message = esc(exs[1])
    kwargs = map(make_kwargs, exs[2:end])
    MPI.Comm_rank(comm) > 0 ? nothing : Expr(:macrocall, Base.CoreLogging.var"@logmsg", location, level, message, kwargs...)
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

