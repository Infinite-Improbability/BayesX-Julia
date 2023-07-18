using SpectralFitting
using Interpolations
using Unitful
using ProgressMeter

"""Mekal model using SpectralFitting.jl framework."""
@xspecmodel :C_mekal struct XS_Mekal{T,F} <: SpectralFitting.AbstractSpectralModel{T,SpectralFitting.Additive}
    "Normalisation"
    K::T
    "Plasma Temperature"
    t::T
    "Hydrogen Density"
    ρ::T
    "Metal Abundances"
    a::T
    "Redshift"
    z::T
    "Switch"
    s::Int
end
function XS_Mekal(;
    K=FitParam(1.0),
    t=FitParam(8.0),
    ρ=FitParam(10.0),
    a=FitParam(1.0),
    z=FitParam(0.1),
    s=0
)
    XS_Mekal{typeof(K),SpectralFitting.FreeParameters{(:K, :t, :ρ)}}(
        K, t, ρ, a, z, s
    )
end
SpectralFitting.register_model_data(XS_Mekal, "mekal1.dat", "mekal2.dat", "mekal3.dat", "mekal4.dat", "mekal5.dat", "mekal6.dat")

"""
    prepare_model_mekal(nHcol, redshift, energy_bins, temperatures; densities, normalisation=1)

Create an interpolated alias to the mekal model over specified parameter ranges.

Interpolation has a significant performance improvement over calling the model directly.
"""
function prepare_model_mekal(
    nHcol,
    redshift,
    energy_bins;
    temperatures=1e-30:0.5:150,
    densities=0:20.0:1000,
    normalisation=1.0
)
    @debug "Preparing MEKAL emission model"

    energy_bins = ustrip.(u"keV", collect(energy_bins))

    # TODO: Figure out normalisation

    # Generate transmission fractions
    @debug "Invoking absorption model"
    absorption_model = PhotoelectricAbsorption(FitParam(nHcol))
    absorption = invokemodel(energy_bins, absorption_model)

    @assert all(isfinite, absorption)

    # Generate source flux
    # TODO: document unit
    @debug "Setting evaluation points"
    # TODO: progress bar
    emission_model(kbT, ρ) = XS_Mekal(K=FitParam(normalisation), t=FitParam(kbT), ρ=FitParam(ρ), z=FitParam(redshift))
    points = [emission_model(t, d) for t in temperatures, d in densities]
    @debug "Invoking MEKAL"
    emission = @showprogress map(
        x -> invokemodel(energy_bins, x),
        points
    )

    @assert all(all.(isfinite, emission))

    # Apply absorption
    @debug "Applying absorption to MEKAL"
    flux = [absorption .* emission[index] for index in eachindex(IndexCartesian(), emission)]

    @assert all(all.(isfinite, flux))

    @debug "Generating interpolation"
    interpol = scale(interpolate(flux, BSpline(Linear())), (temperatures, densities))

    @debug "Emission model generation complete."
    return interpol
end
