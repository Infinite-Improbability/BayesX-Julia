using Unitful, UnitfulAstro

export ConstantModel


"""
    ConstantModel(max_radius, density, temperature)

Constant/uniform density and density model in sphere of given `radius`.

If units are not given radius is assumed to be in kpc, temperature in keV and density in Msun/kpc^3.

This model is intended for testing purposes and cannot accurately replicate a real cluster.
"""
struct ConstantModel{R<:Unitful.Length,T<:Unitful.Temperature,D<:Unitful.Density} <: ClusterModel
    max_radius::R
    density::D
    temperature::T
end
ConstantModel(max_radius::Real, density::Real, temperature::Real)::ConstantModel = ConstantModel(Quantity(max_radius, u"kpc"), Quantity(density, u"Msun/kpc^3"), Quantity(temperature, u"keV"))
# TODO: These unitless fallback constructors annoy me because we make assumptions. Can we find a way to strip them from the priors then dynamically reapply them to the models?

function temperature(cm::ConstantModel, r::Unitful.Length{<:Real})
    r < radius ? cm.temperature : 0.0u"keV"
end

function density(cm::ConstantModel, r::Unitful.Length{<:Real})
    r < radius ? cm.density : 0.0u"g/cm^3"
end
