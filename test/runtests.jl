push!(LOAD_PATH, "src/")
using BayesJ
using MPI
MPI.Init()

include("test_clusters.jl")