using BayesJ
using MPI
using Test
MPI.Init()

default = false
all_tests = false
if isempty(ARGS)
    default = true
elseif "all" in ARGS
    all_tests = true
end

if default || all_tests || "mekal" in ARGS
    include("test_mekal.jl")
end

if default || all_tests || "clusters" in ARGS
    include("test_clusters.jl")
end

if default || all_tests || "background" in ARGS
    include("test_background.jl")
end