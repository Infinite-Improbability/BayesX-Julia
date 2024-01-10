export Model_Constant

"""
    Model_Constant(T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

Constant temperature and density model.

Intended for testing purposes.
"""
function Model_Constant(T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

    temperature(::Unitful.Length{<:Real}) = T
    density(::Unitful.Length{<:Real}) = ρ

    return temperature, density
end

"""
    Model_Constant(T::Unitful.Energy, ρ::Unitful.Density; kwargs...)

Unitless wrapper for [`Model_Constant`](@ref)

Temperature is in keV, density is in g/cm^3
"""
function Model_Constant(T::Real, ρ::Real; kwargs...)
    Model_Constant(
        T * u"keV",
        ρ * u"g/cm^3";
        kwargs...)
end