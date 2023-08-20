# SpectralFitting outputs the current abundance tables when loaded
# But we don't use them so it is misleading
# Redirect output to hide it.
so = stdout
redirect_stdout(devnull)
using SpectralFitting
redirect_stdout(so)

using Interpolations
using Unitful, UnitfulAstro
using ProgressMeter
using LibXSPEC_jll

@derived_dimension SurfaceDensity Unitful.𝐋^-2
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension NumberDensityRate Unitful.𝐋^-3 / Unitful.𝐓

"""Mekal model using SpectralFitting.jl framework. Kept around so we can borrow the model data downloading functions."""
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
    K=FitParam(1.0), # these values were picked arbitarily
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
    call_mekal(energy_range, temperature, nH)

Given a unitful range of energy and unitless temperature (keV) and hydrogen density in the source (cm^-3)
calls MEKAL to calculate the volume emissivity of the source in the source frame.

Due to performance concerns I recommend calling `prepare_model_mekal` instead to generate an interpolation over
the model. The interpolation also implicitly includes absorption, which this function does not.

## Implementation Details

The XSPEC and MEKAL documentation is at times short on detail. Here's what I've figured out, mostly by
inspecting the models in the XSPEC source code. The MEKA model has similar code and more extensive comments that MEKAL
 that are very helpful.

The values calculated by MEKAL are multiplied by a constant of 2.53325e-3. Inspecting the MEKA comments we find the following
equation for this constant which we will term `C`
```math
C = √2 ⋅ h^2 ⋅ α^3 / ((3πm_e)^1.5 ⋅ √(1000e) ⋅ π ⋅ (1pc)^2 ⋅ 10^4)
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

Attempts were made to invoke MEKAL through XSPEC's wrapper functions but I found the modifications they introduced were troublesome to work out,
though some were necessary and are reimplemented in `surface_brightness`. MEKAL's operations are more clearly physically motivated.
"""
function call_mekal(
    energy_range, #keV
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

    if (temperature == 0.0) || (nH == 0.0)
        return fill(0.0u"m^-3/s", n_energy_bins)
    end

    # Scale by volume of and distance to emitting area.
    # MEKAL expects units of 1e50cm^3 / 1pc^2
    # By keeping it at one we match the behaviour of MEKA when
    # it is asked to output volume emissivity rather than a spectrum
    # and so can carry across its coefficents.
    cem = 1.0

    # Abundances of elements w.r.t solar values
    # Using one because it matches BayesX
    # These appear to be hardcoded and thus may not match what XSPEC reports on launch.
    abundances = ones(Cfloat, 15)

    # Initalise output variables
    flux = ones(Cfloat, n_energy_bins)
    ne = 10.0

    # When debugging with GDB it can be helpful to have a breakpoint before entering the
    # Fortran code
    # ccall(:jl_breakpoint, Cvoid, (Any,), flux)

    # We call the FMEKAL function from the XSPEC binary bundled in LibXSPEC_jll
    # The name is mangled to fmekal_ by the compiler and was uncovered
    # using `nm -D libXSFunctions.so` on the shared object, which was found in the
    # Julia artifacts directory.
    # There is a risk this will fail on other platforms as the binaries may be made with
    # different compilers with different mangling conventions, even using the same LibXSPEC_jll version.
    # It may be possible to do some kind of dynamic reference to avoid this?
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

    # We cancel out the 2.53e-3 coefficent applied by MEKAL as it outputs a spectrum in photons/cm^2/s/keV
    # Then we apply the 3.03e-9 coefficent which cancels out the distance-scaled emitting volume `cem` 
    # This converts our results to photons/m^3/s/keV
    # So finally we multiply by the width of our energy bins to get photons/m^3/s/bin
    return 1u"m^(-3)/s/keV" * flux * 3.03103e-9 / 2.53325e-3 .* (max_energy - min_energy)u"keV"
end


"""
    prepare_model_mekal(nHcol, energy_bins, z; temperatures, hydrogen_densities, use_interpolation=true)

Create an interpolated alias to the mekal model with photoelectric absorption over specified parameter ranges.

The interpolation object `(T::Unitful.Energy, nH::Unitful.𝐋^-3) -> volume emissivity` returned gives the volume emissivity per bin of
a source of the specified temperature and hydrogen number density with the effects of absorption from passage through an area
with a hydrogen column density equal to that specified in the `prepare_model_mekal` call.

Applies redshift and time dilation.

Interpolation has a significant performance improvement over calling the model directly.
No redshift is currently applied to energy bins - they should be assumed to be in the source frame.
It may be wise to adjust this but then `surface_brightness` will need adjustment so it doesn't apply the
correction twice.
"""
function prepare_model_mekal(
    nHcol::SurfaceDensity,
    energy_bins::AbstractRange{T},
    z::Real;
    temperatures::AbstractRange{U}=(0:0.05:9.0)u"keV",
    hydrogen_densities::AbstractRange{V}=(0:0.001:1.0)u"cm^-3",
    use_interpolation::Bool=true
) where {T<:Unitful.Energy,U<:Unitful.Energy,V<:NumberDensity}
    @mpidebug "Preparing MEKAL emission model"

    if MPI.Comm_rank(comm) == 0
        @mpidebug "Checking for model data"
        SpectralFitting.download_model_data(PhotoelectricAbsorption, verbose=false, progress=true)
        SpectralFitting.download_model_data(XS_Mekal, verbose=false, progress=true)
    end
    MPI.Barrier(comm)

    # Generate transmission fractions
    @mpidebug "Invoking absorption model"
    absorption_model = PhotoelectricAbsorption(FitParam(ustrip(u"cm^-2", nHcol) / 1e22))
    absorption = invokemodel(ustrip.(u"keV", collect(energy_bins)), absorption_model)

    @assert all(isfinite, absorption)

    # We approximate the absorping hydrogen as at rest w.r.t. observer
    # But the cluster isn't
    # This means we need to apply redshift to the observed energy bins and
    # time dilation to the count rate
    energy_bins = energy_bins * (1 + z) # redshift
    absorption ./= (1 + z) # time dilation

    if !use_interpolation
        function volume_emissivity_direct(
            t::U,
            nH::N
        ) where {U<:Unitful.Energy{Float64},N<:NumberDensity{Float64}}
            let
                energy_bins = energy_bins
                absorption = absorption
                return absorption .* call_mekal(energy_bins, ustrip(u"keV", t), ustrip(u"cm^-3", nH))
            end
        end
        return volume_emissivity_direct
    end

    # Generate source flux
    @mpidebug "Setting evaluation points"
    points = [(ustrip(u"keV", t), ustrip(u"cm^-3", nH)) for t in temperatures, nH in hydrogen_densities]
    @mpidebug "Invoking MEKAL"
    emission = @showprogress 1 "Pregenerating emissions with MEKAL" map(
        x -> replace(ustrip.(u"m^(-3)/s", call_mekal(energy_bins, x...)), NaN => 0),
        points
    )

    @assert all(all.(isfinite, emission))

    # Apply absorption
    @mpidebug "Applying absorption to MEKAL"
    flux = [absorption .* emission[index] for index in eachindex(IndexCartesian(), emission)]

    @assert all(all.(isfinite, flux))

    # The interpolation library doesn't like units so we have to strip and reapply them
    # Also it only wants ranges not vectors for the knots

    @mpidebug "Generating interpolation"
    interpol = Interpolations.scale(interpolate!(flux, BSpline(Linear())), temperatures, hydrogen_densities)

    @mpidebug "Emission model generation complete."

    function volume_emissivity(
        t::U,
        nH::N
    ) where {U<:Unitful.Energy{Float64},N<:NumberDensity{Float64}}
        let
            interpol = interpol
            energy_bins = energy_bins
            absorption = absorption
            try
                return interpol(t, nH) * 1u"m^(-3)/s"
            catch e
                if isa(e, BoundsError)
                    t = uconvert(u"keV", t)
                    nH = uconvert(u"cm^-3", nH)
                    @mpirankeddebug "Exceeded MEKAL interpolation bounds. Calculating the result directly. This is expensive, consider increasing bounds." t nH
                    return absorption .* call_mekal(energy_bins, ustrip(u"keV", t), ustrip(u"cm^-3", nH))
                else
                    throw(e)
                end
            end
        end
    end

    return volume_emissivity
end