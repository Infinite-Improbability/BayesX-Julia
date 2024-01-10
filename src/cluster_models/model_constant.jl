export Model_Constant

"""
    Model_Constant(radius::Unitful.Length, T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

Constant temperature and density model in sphere of given `radius`.

Intended for testing purposes.
"""
function Model_Constant(radius::Unitful.Length, T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

    temperature(r::Unitful.Length{<:Real}) = r < radius ? T : 0.0u"keV"
    density(r::Unitful.Length{<:Real}) = r < radius ? ρ : 0.0u"g/cm^3"

    return temperature, density
end

"""
    Model_Constant(T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

Unitless wrapper for [`Model_Constant`](@ref)

Radius is in Mpc, temperature is in keV, density is in g/cm^3
"""
function Model_Constant(radius::Real, T::Real, ρ::Real; kwargs...)
    Model_Constant(
        radius * 1u"Mpc",
        T * u"keV",
        ρ * u"g/cm^3";
        kwargs...)
end