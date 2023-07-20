using MPI
MPI.Init()

"""Wrappers around the logging macros to ensure they only run for a single process.."""

const comm = MPI.COMM_WORLD

macro mpidebug(exs...)
    return :(MPI.Comm_rank(comm) != 1 ? nothing : @debug $(exs...))
end

macro mpiinfo(exs...)
    return :(MPI.Comm_rank(comm) != 1 ? nothing : @info $(exs...))
end

macro mpiwarn(exs...)
    return :(MPI.Comm_rank(comm) != 1 ? nothing : @warn $(exs...))
end

macro mpierror(exs...)
    return :(MPI.Comm_rank(comm) != 1 ? nothing : @error $(exs...))
end

macro mpirankeddebug(exs...)
    return :(@debug "[Process $(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] $($exs[1] )" $(esc.(exs[2:end])...))
end
