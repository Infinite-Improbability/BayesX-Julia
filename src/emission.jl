using Integrals
using LinearAlgebra: dot
using DimensionfulAngles

include("mekal.jl")

export HydrogenDensity

# Abundances of elements, relative to values from Table 2, Anders & Grevesse (1989) https://doi.org/10.1016/0016-7037(89)90286-X
# H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni
const anders_log_abundance = (12.00, 10.99, 8.56, 8.05, 8.93, 8.09, 6.33, 7.58, 6.47, 7.55, 7.21, 6.56, 6.36, 7.67, 6.25)
const ander_Ni_per_NH = 10 .^ (anders_log_abundance .- 12)
const nucleon_total = (1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0, 28.0, 32.0, 40.0, 40.0, 56.0, 59.0)

"""
    ObservationError(likelihood)

This wraps a fallback likelihood value so it can be passed up the chain on an invalid
prior combination.
"""
struct ObservationError <: Exception
    likelihood::Float64
end

"""
    surface_brightness(projected_radius, temperature, density, z, limit, model, pixel_edge_angle)

Calculate the observed surface_brightness at some projected radius on the sky.

The temperature and density functions should take the distance from origin (cluster centre) and return the appropriate value, with units.

The model should take these functions as input and return a vector of volume emissivities, where the vector element correspond with
the energy bins used to generate the mode. It should include the effects of redshift on the bins and time dilation on the count rate.

The limit controls limits on LOS radius during integration. Testing suggests that Inf (default) is faster than finite values >1.
However finite values may avoid problems with models that are badly behaved at large radii.

The pixel edge angle is the angular size of a single pixel.

The result returned is the expected count rate per unit observing area, as a vector corresponding to the energy bins as used in the model.
"""
function surface_brightness(
    projected_radius::Unitful.Length,
    temperature,
    density,
    z::Float64,
    limit::Unitful.Length,
    model!::Function,
    pixel_edge_angle::DimensionfulAngles.Angle,
    flux::Vector{Float32}
)::Vector{Quantity{Float64,Unitful.𝐋^(-2) / Unitful.𝐓}}
    @argcheck limit > 0u"Mpc"

    lim = ustrip(Float64, u"m", limit)
    pr = ustrip(Float64, u"m", projected_radius)

    flux .= 0.0f0

    function integrand(y::Vector{Float32}, l::AbstractFloat, params::Tuple{Float64,Any,Any})
        r = Quantity(hypot(params[1], l), u"m")

        # Testing shows that swapping to explicitly Mpc^-3 s^-1 makes ~1e-14 % difference to final counts
        # Result is in m^-3/s
        model!(y, params[2](r), params[3](r))
    end

    # Only integrate from 0 to limit because it is faster and equal to 1/2 integral from -limit to limit
    ifunc = IntegralFunction(integrand, flux)
    problem = IntegralProblem(ifunc, (0.0, lim), (pr, temperature, density))
    sol = solve(problem, HCubatureJL(); reltol=1e-2)
    u = sol.u * 1u"m^-2/s"
    if all(isfinite, sol.u) == false
        @mpirankedwarn "Integration returned non-finite values. Returning fallback likelihood." projected_radius sol.u
        throw(ObservationError(-1e100 * (length(sol.u) - count(isfinite.(sol.u)))))
    elseif all(iszero, sol.u)
        @mpirankeddebug "Integration found point without emission" projected_radius temperature(0u"kpc") temperature(1u"kpc") temperature(10u"kpc") temperature(100u"kpc") density(0u"kpc") density(1u"kpc") density(10u"kpc") density(100u"kpc")
    end
    return 2 * u * pixel_edge_angle^2 / Quantity(4π, u"srᵃ") / (1 + z)^2

    # sol is volume emissivity per face area of column
    # because of how we defined the limits we have to double it
    # [photons/s/m^2]
    # u = (2 * sol.u)

    # add in surface area of column end
    # the pixel is square with angular size θ = pixel_edge_angle
    # this converts to the physical length l = dₐ * θ where dₐ is the angular diameter distance
    # The true area is thus (observed angular area * angular diamemter distance)^2
    # [photons/s]
    # u = (2 * sol.u) * (dₐ * θ)²

    # time dilation and redshift are already factored into model
    # because redshift needs to be applied directly to the energy bins
    # and time dilation felt more connected to that than the spatial expansion at play here

    # get emission per solid angle, assuming uniformly distributed over a sphere
    # [photons/s/sr]
    # u = (2 * sol.u) * (dₐ * θ)² / 4π

    # The radius of the sphere is dₘ which is the transverse comoving distance
    # dₐ = dₘ / (1+z)
    # So we can cancel out the dₐ
    # [photons/s/sr/m^2]
    # u = (2 * sol.u) * (dₐ * θ)² / 4π / dₘ^2
    # u = (2 * sol.u) * (dₘ / (1+z) * θ)² / 4π / dₘ^2
    # u = (2 * sol.u) * θ² / 4π / (1+z)²

end

"""
    apply_response_function(counts_per_bin::Vector, response::Matrix, exposure_time::Unitful.Time)::Vector{Float64}

Applies the response function ``RSP(PI, E) = RMF(PI, E) ∘ ARF(E)`` to all energy bins and returns adjusted counts per bin.
```math
\\begin{aligned}
C(PI) &= T \\int RMF(PI, E) ⋅ ARF(E) ⋅ S(E) ⋅ dE \\\\
&≈ T \\sum_{j} R_{ij} A_{j} S_{j}
\\end{aligned}
```
with `C(PI)` is the observed counts in a detector channel `PI`, `T` is the observation time,
`ARF(E)` is the effective area of the instrument and `RMF(E, PI)` [we use `RMF(PI, E)`] is the unitless response matrix.
`S(E)` is the source flux, in [photon count]/[time]/[area]/[energy].
(Handbook of X-ray Astronomy Chapter 2011, Chapter 5 pg 87, by K. Arnaud, R. Smith and A. Siemiginowska)

This function takes the combined RMF and ARF as the response function. This is to recalculating it on every call.
Some people format the RMF as RMF(E, PI), others as RMF(PI, E). This latter convention is used by CIAO and by us.
"""
function apply_response_function(counts_per_bin::Vector{T}, response::Matrix{T}, exposure_time::T)::Vector{T} where {T<:AbstractFloat}
    @argcheck length(counts_per_bin) == size(response)[2]

    for i in eachindex(counts_per_bin)
        counts_per_bin[i] = counts_per_bin[i] * exposure_time
    end

    # All we have to do is matrix multiplication
    mult = response * counts_per_bin
    return mult
end

struct HydrogenDensity{T<:AbstractVector{<:Real},U<:Unitful.Mass}
    relative_abundances::T
    gas_mass_per_hydrogen::U
end
function HydrogenDensity(relative_abundances::AbstractVector{<:Real})
    # Apply our weights to the Anders and Grevesse abundances
    Ni_per_NH = ander_Ni_per_NH .* relative_abundances

    # Precalculate some values
    total_nucelons_per_hydrogen = dot(Ni_per_NH, nucleon_total)
    gas_mass_per_hydrogen = m_p * total_nucelons_per_hydrogen

    HydrogenDensity(relative_abundances, gas_mass_per_hydrogen)
end
function (gm::HydrogenDensity)(ρ::Unitful.Density)
    return ρ / gm.gas_mass_per_hydrogen
end

function hydrogen_number_density(gas_density::Unitful.Density, relative_abundance=ones(15))
    # This function has been replaced with the form above for performance reasons - precalculated values.
    # Abundances of elements, relative to values from Table 2, Anders & Grevesse (1989) https://doi.org/10.1016/0016-7037(89)90286-X
    # H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni
    abundance = (12.00, 10.99, 8.56, 8.05, 8.93, 8.09, 6.33, 7.58, 6.47, 7.55, 7.21, 6.56, 6.36, 7.67, 6.25)
    Ni_per_NH = 10 .^ (abundance .* relative_abundance .- 12)
    nucleon_total = (1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0, 28.0, 32.0, 40.0, 40.0, 56.0, 59.0)
    return gas_density / (m_p * dot(Ni_per_NH, nucleon_total))
end
# @deprecate hydrogen_number_density(gas_density::Unitful.Density, relative_abundance::AbstractVector{<:Real}) HydrogenDensity(relative_abundance)(gas_density::Unitful.Density)

"""
    pixel_offset(ij, array_centre_indices, centre_offset_pixels)

Calculate the distance from the centre of the cluster to the pixel at indices ij.

Result in pixels.
"""
function pixel_offset(ij, array_centre_indices, centre_offset_pixels)::Float64
    offset = Tuple(ij) .- array_centre_indices .- centre_offset_pixels
    return hypot(offset...)
end

"""
    make_observation(temperature, density, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre, centre_radius, mask=nothing, limit=10u"Mpc")

Generate an image of the cluster given functions for the radial profile of gas temperature and electron density and assorted observational parameters.

- `shape` describes the size of the source number as a number of spatial bins in each dimension.-
- The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds. This area is assumed to be square with the edge angle giving the side length.
- The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
- The response function includes both the RMF and ARF, as described in [`apply_response_function`](@ref).
- The limit is passed through to [`surface_brightness`](@ref)
"""
function make_observation(
    temperature,
    density,
    z::Real,
    shape::NTuple{2,<:Integer},
    pixel_edge_angle::A,
    emission_model,
    exposure_time::T,
    response_function::AbstractArray{<:Unitful.Area,2},
    centre::NTuple{2,<:DimensionfulAngles.Angle},
    centre_radius::Real;
    mask::Union{Matrix{Bool},Nothing}=nothing,
    limit::Unitful.Length=Quantity(10, u"Mpc")
)::Array{Union{Float64,Missing},3} where {A<:DimensionfulAngles.Angle,T<:Unitful.Time}
    # pixel size
    pixel_edge_length = angle_to_length(pixel_edge_angle, z)

    # Centre of the data array
    array_centre_pixels = (shape .+ 1) ./ 2

    # Offset of cluster centre from data array centre
    centre_offset = angle_to_length.(centre, z)
    centre_offset_pixels = centre_offset ./ pixel_edge_length
    core_position = array_centre_pixels .+ centre_offset_pixels

    if any(centre_radius .> array_centre_pixels)
        @mpirankedwarn "Centre exclusion radius greater than observed radius in at least one direction" centre_radius array_centre_pixels
    end

    brightness_radii = range(
        centre_radius,
        stop=hypot((core_position .+ 2)...),
        step=2.0
    )
    @mpirankeddebug "Creating brightness interpolation" length(brightness_radii)

    flux = Vector{Float32}(undef, size(response_function, 2))
    brightness_line = Vector{Vector{Float32}}(undef, length(brightness_radii))
    for i in eachindex(brightness_line)
        brightness_line[i] = ustrip.(
            Float64,
            u"m^(-2)/s",
            surface_brightness(
                brightness_radii[i] * pixel_edge_length,
                temperature,
                density,
                z,
                limit,
                emission_model,
                pixel_edge_angle,
                flux
            )
        )
    end

    brightness_interpolation = linear_interpolation(brightness_radii, brightness_line, extrapolation_bc=Throw())

    @mpirankeddebug "Calculating counts"
    @mpirankeddebug "Preparing response function"
    resp::Matrix{Float64} = ustrip.(Float64, u"m^2", response_function)
    exp_time::Float64 = ustrip(Float64, u"s", exposure_time)
    @mpirankeddebug "Creating counts array"
    counts::Array{Union{Float64,Missing}} = fill(missing, size(response_function, 1), shape...)

    @mpirankeddebug "Checking mask"
    if isnothing(mask)
        mask = zeros(Bool, shape...)
    end

    possible_indices = CartesianIndices(ceil.(Int, shape))
    possible_radii = unique((pixel_offset(ij, array_centre_pixels, centre_offset_pixels) for ij in possible_indices))
    possible_radii = filter!(r -> r >= centre_radius, possible_radii)

    emissions = Dict{Float64,Vector{Float64}}(r => apply_response_function(brightness_interpolation(r), resp, exp_time) for r in possible_radii)

    @mpirankeddebug "Iterating over pixels"
    for ij in CartesianIndices(shape)
        if !mask[ij[1], ij[2]]
            radius = pixel_offset(ij, array_centre_pixels, centre_offset_pixels)
            if radius >= centre_radius
                counts[:, ij[1], ij[2]] = emissions[radius]
            end
        end
    end

    # @mpirankeddebug "Checking counts"

    # if all(ismissing, counts)
    #     @mpierror "All counts generated by make_observation are missing. This is likely due to masking problems"
    # elseif all(iszero, counts) !== false # can return true, false or missing (if all non-missing values are true)
    #     # this should be unreachable with the earlier checks
    #     @mpidebug "All counts generated by make_observation are zero. Somehow."
    # end

    return counts
end
"""
    make_observation(temperature, density, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre, centre_radius, mask=nothing)

Unitless wrapper for [`make_observation`](@ref)
"""
function make_observation(
    temperature,
    density,
    z::Real,
    shape::NTuple{2,<:Integer},
    pixel_edge_angle::A,
    emission_model,
    exposure_time::T,
    response_function,
    centre::NTuple{2,Real},
    centre_radius;
    mask::Union{Matrix{Bool},Nothing}=nothing,
    limit::Unitful.Length=Quantity(Inf, u"Mpc")
)::Array{Union{Float64,Missing},3} where {A<:DimensionfulAngles.Angle,T<:Unitful.Time}
    @mpidebug "Called make_observation wrapper"
    make_observation(
        temperature,
        density,
        z,
        shape,
        pixel_edge_angle,
        emission_model,
        exposure_time,
        response_function,
        centre .* 1u"arcsecondᵃ",
        centre_radius;
        mask=mask,
        limit=limit
    )
end
