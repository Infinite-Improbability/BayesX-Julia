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

"""
    call_mekal(energy_range, temperature, nH)

Given a unitful range of energy and unitless temperature (keV) and hydrogen density in the source (cm^-3)
calls MEKAL to calculate the volume emissivity of the source.

The XSPEC and MEKAL documentation is at times short on detail. Here's what I've figured out.

The MEKA model has more extensive comments that MEKAL that are very helpful. The values calculated
by MEKAL are multiplied by a constant of 2.53325e-3. Inspecting the MEKA comments we find the following
equation for this constant which we will term `C`
```math
C = \\frac{\\sqrt{2} h^2 α^3}{(3πm_e)^{1.5} \\sqrt{1000e} π (1pc)^2 10^4}
```
and the instruction to use SI units. `h` appears to be Planck's constant (unreduced) and `α` the fine-structure constant.
Calculating this gives C=2.533242122484874e-59 kg^1/2 m^2 A^-1/2 s^-5/2.

MEKAL adjusts emissions by what it terms `cem` which is the emitting volume divided by the square of distance to the volume.
This follows naturally from the relationship between surface emissivity `σ` and volume emissivity `ϵ`.
```math
    σ = ϵV / (4πD)^2
```
MEKAL uses units of 1e50cm^3 / 1pc^2 for `cem`. If we multiply by 10^50cm^3 we get 2.533242122484874e-13 J^2 cm s^2 C^-1/2 kg^-3/2.
As dividing by 1pc^2 would change more than just the magnitude of the value I assume it is folded in as the 1pc^2 in the original equation.
Furthermore I assume that remaining difference rests in the details of units and implementation.

MEKAL states it outputs in phot/cm^2/s/keV. MEKA states the same but has an optional switch that instead outputs in phot/m^3/s/keV.
It does this by replacing `C` with `D=3.03103E-9` which it states is derived so that
```math
C = D * 10^40 / (4π * 1pc^2)
```
I have verified this gives results matching that of doing the conversion manually (by dividing out 1e50cm^3 and multiplying by 1pc^2).
"""
function call_mekal(
    energy_range,
    # emitting_volume,
    # distance,
    temperature, # keV
    nH, # cm^-3
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

    return 1u"m^(-3)/s/keV" * flux * 3.03103e-9 / 2.53325e-3 .* (max_energy - min_energy)u"keV"
end

# display(call_mekal((0.3:0.1:7.0)u"keV", 10, 1e-3))

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