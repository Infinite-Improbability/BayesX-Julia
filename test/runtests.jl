push!(LOAD_PATH, "src/")
using BayesJ
using MPI
MPI.Init()

# ENV["JULIA_DEBUG"] = "BayesJ"

include("test_clusters.jl")