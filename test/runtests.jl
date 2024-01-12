using BayesJ
using MPI
using Test
MPI.Init()

if isempty(ARGS) || "all" in ARGS
    all_tests = true
else
    all_tests = false
end

if all_tests || "clusters" in ARGS
    include("test_clusters.jl")
end

if all_tests || "mekal" in ARGS
    include("test_mekal.jl")
end

if all_tests || "internal" in ARGS
    include("test_internal_consistency.jl")
end