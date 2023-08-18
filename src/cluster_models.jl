using ArgCheck
using Unitful, UnitfulAstro, DimensionfulAngles
using PhysicalConstants.CODATA2018: G
using Integrals

include("params.jl")
include("emission.jl")

"""
    p_crit(z)

Calculate the critical density at some redshift `z`.
"""
ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)


include("cluster_models/model_nfw.jl")
include("cluster_models/model_einasto.jl")
include("cluster_models/model_vikhlinin2006.jl")
