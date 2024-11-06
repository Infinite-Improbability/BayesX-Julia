using BenchmarkTools
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using LinearAlgebra: I
using Profile

redshift = 0.1
nfw = NFWModel(3e14u"Msun", 0.13, 4.0, 1.0620, 5.4807, 0.3292, 1.156, z=redshift)
energy_range = range(0.1u"keV", 13u"keV", length=100)
mekal = prepare_model_mekal(
    0.022e22u"cm^-2",
    energy_range,
    redshift
)
response_function = Matrix{Float64}(I, length(energy_range) - 1, length(energy_range) - 1) * 1u"cm^2"

function make_obs_wrapper()
    BayesJ.make_observation(
        nfw,
        redshift,
        (30, 30),
        0.492u"arcsecondᵃ",
        mekal,
        3000u"s",
        response_function,
        (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        1,
    )
end

make_obs_wrapper(); # force compilation

display(@benchmark make_obs_wrapper())
@profview make_obs_wrapper();
