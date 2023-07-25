using SpectralFitting
using Interpolations
using Unitful, UnitfulAstro
using ProgressMeter
using LibXSPEC_jll

libXSFunctions

include("mpi.jl")

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
    prepare_model_mekal(nHcol, redshift, energy_bins; temperatures, densities, normalisation=1)

Create an interpolated alias to the mekal model over specified parameter ranges.

Interpolation has a significant performance improvement over calling the model directly.
"""
function prepare_model_mekal(
    nHcol,
    redshift,
    energy_bins;
    temperatures=1e-30:5.0:500.0,
    densities=0:5.0:200,
    normalisation=1.0
)
    @mpidebug "Preparing MEKAL emission model"

    energy_bins = ustrip.(u"keV", collect(energy_bins))

    # TODO: Figure out normalisation

    # Generate transmission fractions
    @mpidebug "Invoking absorption model"
    absorption_model = PhotoelectricAbsorption(FitParam(nHcol))
    absorption = invokemodel(energy_bins, absorption_model)

    @assert all(isfinite, absorption)

    # Generate source flux
    # TODO: document unit
    @mpidebug "Setting evaluation points"
    # TODO: progress bar
    emission_model(kbT, ρ) = XS_Mekal(K=FitParam(normalisation), t=FitParam(kbT), ρ=FitParam(ρ), z=FitParam(redshift))
    points = [emission_model(t, d) for t in temperatures, d in densities]
    @mpidebug "Invoking MEKAL"
    emission = @showprogress 1 "Pregenerating emissions with MEKAL" map(
        x -> invokemodel(energy_bins, x),
        points
    )

    @assert all(all.(isfinite, emission))

    # Apply absorption
    @mpidebug "Applying absorption to MEKAL"
    flux = [absorption .* emission[index] for index in eachindex(IndexCartesian(), emission)]

    @assert all(all.(isfinite, flux))

    @mpidebug "Generating interpolation"
    interpol = scale(interpolate(flux, BSpline(Linear())), (temperatures, densities))

    @mpidebug "Emission model generation complete."
    return interpol
end


function call_mekal(
    energy_range,
    # emitting_volume,
    # distance,
    temperature,
    nH,
)
    # Convert energy range into format expected by mekal
    n_energy_bins = length(energy_range) - 1
    min_energy = Vector{Cfloat}(undef, n_energy_bins)
    max_energy = Vector{Cfloat}(undef, n_energy_bins)
    for i in 1:n_energy_bins
        min_energy[i] = ustrip(Cfloat, u"keV", energy_range[i])
        max_energy[i] = ustrip(Cfloat, u"keV", energy_range[i+1])
    end

    # Scale by volume of and distance to emitting area.
    # cem_units = uconvert(u"m", 1e50u"cm^3" / (1u"pc")^2)
    # cem = uconvert(u"m", emitting_volume / distance^2) / cem_units
    cem = 1.0

    # Ensure input quantities are in the right units
    # nH = ustrip(Cfloat, u"cm^-3", nH)
    # temperature = ustrip(Cfloat, u"keV", temperature)

    # Abundances of elements w.r.t solar values
    # Using one because it matches BayesX
    abundances = ones(Cfloat, 15)

    # Initalise output variables
    flux = ones(Cfloat, n_energy_bins)
    ne = 10.0

    # ccall(:jl_breakpoint, Cvoid, (Any,), flux)

    @ccall "libXSFunctions".fmekal_(
        min_energy::Ptr{Cfloat},
        max_energy::Ptr{Cfloat},
        flux::Ptr{Cfloat},
        n_energy_bins::Ref{Cint},
        cem::Ref{Cfloat},
        nH::Ref{Cfloat},
        temperature::Ref{Cfloat},
        abundances::Ptr{Cfloat},
        ne::Ref{Cfloat},
    )::Cvoid

    return 1u"m^(-3)/s/keV" * 3.03103e-9 / 2.53325e-3 * flux .* (max_energy - min_energy)u"keV"

    # return 1e-4u"m^(-3)/s/keV" * flux * 3.03103e-9 / 2.53325e-3 .* (max_energy - min_energy)u"keV"
end

"""
    prepare_model_mekal2(nHcol, energy_bins; temperatures, densities)

Create an interpolated alias to the mekal model over specified parameter ranges.

Interpolation has a significant performance improvement over calling the model directly.
This version calls MEKAL without using XSPEC's wrappers.
"""
function prepare_model_mekal2(
    nHcol,
    energy_bins;
    temperatures=1e-30:5.0:500.0, #keV
    densities=0:5.0:200 # cm^-3
)
    @mpidebug "Preparing MEKAL emission model"

    # TODO: Do I need to redshift the energy bins?

    # Generate transmission fractions
    @mpidebug "Invoking absorption model"
    absorption_model = PhotoelectricAbsorption(FitParam(nHcol))
    absorption = invokemodel(ustrip.(u"keV", collect(energy_bins)), absorption_model)

    @assert all(isfinite, absorption)

    # Generate source flux
    @mpidebug "Setting evaluation points"
    # TODO: progress bar
    points = [(t, d) for t in temperatures, d in densities]
    @mpidebug "Invoking MEKAL"
    emission = @showprogress 1 "Pregenerating emissions with MEKAL" map(
        x -> ustrip.(u"m^(-3)/s", call_mekal(energy_bins, x...)),
        points
    )

    @assert all(all.(isfinite, emission))

    # Apply absorption
    @mpidebug "Applying absorption to MEKAL"
    flux = [absorption .* emission[index] for index in eachindex(IndexCartesian(), emission)]

    @assert all(all.(isfinite, flux))

    @mpidebug "Generating interpolation"
    interpol = scale(interpolate(flux, BSpline(Linear())), temperatures, densities)

    @mpidebug "Emission model generation complete."
    f(t, d) = interpol(ustrip(u"keV", t), ustrip(u"cm^-3", d)) * 1u"m^(-3)/s"
    return f
end