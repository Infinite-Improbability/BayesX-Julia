push!(LOAD_PATH, "src/")
using BayesJ
using MPI
MPI.Init()
using Test

# ENV["JULIA_DEBUG"] = "BayesJ"

include("test_clusters.jl")