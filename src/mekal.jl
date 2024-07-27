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

export prepare_model_mekal, prepare_model_mekal_interpolation

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


const mekal_factor = 3.03103f-9 / 2.53325f-3

"""
    call_mekal(n_energy_bins, abundances, min_energy, max_energy, bin_sizes, temperature, nH)

Given a unitful range of energy and unitless temperature (keV) and hydrogen density in the source (cm^-3)
calls MEKAL to calculate the volume emissivity of the source in the source frame.

Result is in photons/m^3/s/bin.

Metal abundances are set by `abundances`. They include H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni and are specified with respect to solar values.
MEKAL uses the solar values from [andersAbundancesElementsMeteoritic1989](@cite).

Due to performance concerns I recommend calling `prepare_model_mekal` instead to generate an interpolation over
the model. The functions returned by that wrapper also implicitly include absorption, which this function does not.

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
    flux::Vector{Cfloat},
    abundances::Vector{Cfloat},
    n_energy_bins::Integer,
    min_energy::V, #keV
    max_energy::V, #keV
    bin_sizes::V,
    temperature::Cfloat, # keV
    nH::Cfloat, # cm^-3
) where {V<:AbstractVector{Cfloat}}
    if (temperature == 0.0) || (nH == 0.0)
        for i in 1:n_energy_bins
            @inbounds flux[i] = 0.0f0
        end
        return
    end

    # cem: Scale by volume of and distance to emitting area.
    # MEKAL expects units of 1e50cm^3 / 1pc^2
    # By keeping it at one we match the behaviour of MEKA when
    # it is asked to output volume emissivity rather than a spectrum
    # and so can carry across its coefficents.

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
        1.0f0::Ref{Cfloat}, # cem
        nH::Ref{Cfloat},
        temperature::Ref{Cfloat},
        abundances::Ptr{Cfloat},
        10.0f0::Ref{Cfloat}, # ne
    )::Cvoid

    # We cancel out the 2.53e-3 coefficent applied by MEKAL as it outputs a spectrum in photons/cm^2/s/keV
    # Then we apply the 3.03e-9 coefficent which cancels out the distance-scaled emitting volume `cem` 
    # This converts our results to photons/m^3/s/keV
    # So finally we multiply by the width of our energy bins to get photons/m^3/s/bin
    @simd for i in 1:n_energy_bins
        @inbounds flux[i] = flux[i] * bin_sizes[i] * mekal_factor
    end
end
function call_mekal(
    flux::AbstractVector{F},
    abundances::AbstractVector{A},
    energy_range::AbstractVector{E},
    temperature::Real, # keV
    nH::Real, # cm^-3
) where {F<:AbstractFloat,A<:AbstractFloat,E<:Unitful.Energy}

    # Convert energy range into format expected by mekal
    n_energy_bins = length(energy_range) - 1
    min_energy = ustrip.(Cfloat, u"keV", energy_range[1:end-1])
    max_energy = ustrip.(Cfloat, u"keV", energy_range[2:end])
    bin_sizes = max_energy - min_energy

    call_mekal(
        convert(Vector{Cfloat}, flux),
        convert(Vector{Cfloat}, abundances),
        n_energy_bins,
        min_energy,
        max_energy,
        bin_sizes,
        convert(Cfloat, temperature), convert(Cfloat, nH)
    )
end

"""
    download_spectral_fitting_data(; verbose=false, progress=true, kwargs...)

Download SpectralFitting.jl model data.

Keyword arguments passed through to `SpectralFitting.download_model_data`.
"""
function download_spectral_fitting_data(; verbose=false, progress=true, kwargs...)
    SpectralFitting.download_model_data(PhotoelectricAbsorption, verbose=verbose, progress=progress, kwargs...)
    SpectralFitting.download_model_data(XS_Mekal, verbose=verbose, progress=progress, kwargs...)
end


"""
    prepare_model_mekal(nHcol, energy_bins, z, [abundances])

Generate a function to calculate emission with MEKAL .

The returned function has signature
```
f!(flux::Vector{Float32}, T::Unitful.Energy, nH::Unitful.𝐋^-3)::Vector{Float32}
```

It works inplace on `flux`. After execution `flux` contains the volume emissivity in photons per cubic metre per second per bin of 
a source of the specified temperature and hydrogen number density. Photon redshift and time dilation effects based on the redshift `z`. 
It also applies the effects of absorption from passage through an area with a hydrogen column density `nHcol`.

Metal abundances are set by `abundances`. They include H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni and are specified with respect to solar values.
MEKAL uses the solar values from [andersAbundancesElementsMeteoritic1989](@cite).

"""
function prepare_model_mekal(
    nHcol::SurfaceDensity,
    energy_bins::AbstractVector{T},
    z::Real,
    abundances::AbstractVector{A}=ones(15),
) where {A<:Real,T<:Unitful.Energy}
    @mpidebug "Preparing MEKAL emission model"

    @argcheck length(abundances) == 15
    @argcheck abundances[1] == 1

    if MPI.Comm_rank(comm) == 0
        @mpidebug "Checking for model data"
        download_spectral_fitting_data()
        @mpidebug "Model data ready"
    end
    MPI.Barrier(comm)

    # Generate transmission fractions
    @mpidebug "Invoking absorption model"
    absorption_model = PhotoelectricAbsorption(FitParam(ustrip(u"cm^-2", nHcol) / 1e22))
    absorption = invokemodel(ustrip.(u"keV", energy_bins), absorption_model)

    @assert all(isfinite, absorption)
    @assert all(absorption .<= 1.0)
    @assert all(absorption .>= 0.0)

    # We approximate the absorping hydrogen as at rest w.r.t. observer
    # But the cluster isn't
    # This means we need to apply redshift to the observed energy bins and
    # time dilation to the count rate

    # the original bins are at the observer and E_obs = E_emit / (1+z)
    # to get the bins at the source for mekal we multiply by 1+z
    energy_bins = energy_bins * (1 + z) # redshift

    # absorption is applied to the source counts
    absorption ./= (1 + z) # time dilation

    # Convert energy range into format expected by mekal
    n_energy_bins = length(energy_bins) - 1
    min_energy = ustrip.(Cfloat, u"keV", energy_bins[1:end-1])
    max_energy = ustrip.(Cfloat, u"keV", energy_bins[2:end])
    bin_sizes = max_energy - min_energy

    nH = HydrogenDensity(abundances)

    # Use the values relative to Anders and Grevesse
    abundances_float = convert(Vector{Cfloat}, abundances)

    # TODO: Rewrite inner function to use struct, like we've done with hydrogen_number_density
    # That should retain performance without needing to use let so much.

    function volume_emissivity!(
        flux::Vector{Cfloat},
        t::U,
        ρ::N
    ) where {U<:Unitful.Energy,N<:Unitful.Density}
        let n_energy_bins = n_energy_bins, min_energy = min_energy, max_energy = max_energy, bin_sizes = bin_sizes, absorption = absorption, abundances_float = abundances_float, nH = nH
            call_mekal(
                flux,
                abundances_float,
                n_energy_bins,
                min_energy,
                max_energy,
                bin_sizes,
                ustrip(Cfloat, u"keV", t),
                ustrip(Cfloat, u"cm^-3", nH(ρ))
            )
            @simd for i in 1:n_energy_bins
                @inbounds flux[i] = flux[i] * absorption[i]
            end
        end
    end
    return volume_emissivity!


end

"""
prepare_model_mekal_interpolation(nHcol, energy_bins, z, [abundances, temperatures, hydrogen_densities])

Generate an interpolation function to calculate emission with MEKAL.

This wraps [`prepare_model_mekal`](@ref) to generate an interpolation function over the model. 
Interpolation points are determined from `temperatures` and `hydrogen_densities`.

Interpolation has a significant performance improvement over calling the model directly.
No redshift is currently applied to energy bins - they should be assumed to be in the source frame.
It may be wise to adjust this but then `surface_brightness` will need adjustment so it doesn't apply the
correction twice.
"""
function prepare_model_mekal_interpolation(
    nHcol::SurfaceDensity,
    energy_bins::AbstractVector{T},
    z::Real,
    abundances::AbstractVector{A}=ones(15);
    temperatures::AbstractRange{U}=range(0u"keV", 9u"keV", 1000),
    gas_densities::AbstractRange{V}=range(0u"g/cm^3", m_p * 1u"cm^-3", 1000),
) where {A<:Real,T<:Unitful.Energy,U<:Unitful.Energy,V<:Unitful.Density}

    # Get wrapper around standard MEKAL call
    base_model = prepare_model_mekal(nHcol, energy_bins, z, abundances)

    # setup output array
    flux::Vector{Float32} = zeros(Float32, length(energy_bins) - 1)
    function calc_flux(T, nH)
        base_model(flux, T, nH)
        return flux
    end

    # Generate source flux
    total_points = length(temperatures) * length(gas_densities)
    @mpidebug "Setting MEKAL evaluation points" size(temperatures) size(gas_densities) total_points
    points = [(t, nH) for t in temperatures, nH in gas_densities]
    @mpidebug "Invoking MEKAL"
    emission = @showprogress 1 "Pregenerating emissions with MEKAL" map(
        x -> calc_flux(x[1], x[2]),
        points
    )
    @assert all(all.(isfinite, emission))

    @mpidebug "Generating interpolation"
    interpol = Interpolations.scale(interpolate!(emission, BSpline(Linear())), temperatures, gas_densities)
    @mpidebug "Emission interpolation generation complete."

    function volume_emissivity!(
        flux::Vector{Cfloat},
        t::U,
        ρ::N
    ) where {U<:Unitful.Energy,N<:Unitful.Density}
        let interpol = interpol, base_model = base_model
            try
                # return interpolation if possible
                flux .= interpol(t, ρ)
            catch e
                # otherwise fall back to direct call
                if isa(e, BoundsError)
                    base_model(flux, t, ρ)
                else
                    rethrow(e)
                end
            end
        end
    end

    return volume_emissivity!
end