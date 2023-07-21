using MPI
MPI.Init()

"""Wrappers around the logging macros to ensure they only run for a single process.."""

const comm = MPI.COMM_WORLD

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

macro mpirankedinfo(msg, exs...)
    return :(@info "[Process $(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] $($msg)" $(esc.(exs)...))
end


macro mpirankeddebug(msg, exs...)
    return :(@debug "[Process $(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] $($msg)" $(esc.(exs)...))
end
