using Integrals
using LinearAlgebra: dot
using DimensionfulAngles

include("mekal.jl")

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

The limit controls limits on LOS radius during integration. Theoretically it should be ±infinity, but we may approximate it as a finite value.

The pixel edge angle is the angular size of a single pixel.

The result returned is the expected count rate per unit observing area, as a vector corresponding to the energy bins as used in the model.
"""
function surface_brightness(
    projected_radius::Unitful.Length,
    temperature::Function,
    density::Function,
    z::Float64,
    limit::Unitful.Length,
    model,
    pixel_edge_angle::DimensionfulAngles.Angle,
)::Vector{Quantity{Float64,Unitful.𝐋^(-2) / Unitful.𝐓}}
    @argcheck limit > 0u"Mpc"

    lim = ustrip(Float64, u"Mpc", limit)
    pr = ustrip(u"Mpc", projected_radius)

    nout = length(model(1.0u"keV", 0.1u"cm^-3"))

    function integrand(l, params)
        s, temp, density = params
        r = hypot(s, l) * 1u"Mpc"

        # Testing shows that swapping to explicitly Mpc^-3 s^-1 makes ~1e-14 % difference to final counts
        f = model(temp(r), hydrogen_number_density(density(r)))

        return ustrip.(u"Mpc^-3/s", f)
    end

    # Only integrate from 0 to limit because it is faster and equal to 1/2 integral from -limit to limit
    problem = IntegralProblem(integrand, 0.0, lim, (pr, temperature, density); nout=nout)
    sol = solve(problem, HCubatureJL(); reltol=1e-3, abstol=1.0)
    # @assert all(isfinite, sol.u)
    u = sol.u
    if all(isfinite, sol.u) == false
        @mpirankedwarn "Integration returned non-finite values. Returning fallback likelihood."
        raise(ObservationError(-1e100 * (length(sol.u) - count(isfinite(sol.u)))))
    end
    return 2 * u * 1u"Mpc^-2/s" / (Quantity(4π, u"srᵃ") * (1 + z)^2) * pixel_edge_angle^2

    # sol is volume emissivity per face area of column
    # because of how we defined the limits we have to double it
    # [photons/s/m^2]
    # u::Vector{Quantity} = 2 * sol.u

    # add in surface area of column end
    # The true anglular area is (1+z)^2 * observed but using the angular diameter distance should avoid that problem.
    # [photons/s]
    # u *= (dₐ * pixel_edge_angle)^2

    # time dilation and redshift are already factored into model
    # because redshift needs to be applied directly to the energy bins
    # and time dilation felt more connected to that than the spatial expansion at play here

    # get emission per solid angle, assuming uniformly distributed over a sphere
    # u /= 4π u"sr"

    # but this is the solid angle at the source
    # for the observer θ₀ = θₛ / (1 + z)
    # u /= (1+z)^2

    # convert from solid angle to area using angular diameter distance
    # u /= (4π * dₐ^2)
    # notice dₐ cancels out

end

"""
    apply_response_function(counts_per_bin::Vector, response::Matrix, exposure_time::Unitful.Time)::Vector{Float64}

Applies the response function ``RSP(PI, E) = RMF(PI, E) ∘ ARF(E)`` to all energy bins and returns adjusted counts per bin.
```math
\\begin{aligned}
C(PI) &= T \\int RMF(PI, E) ⋅ ARF(E) ⋅ S(E) ⋅ dE \\\\
&≈ T \\sum_{j} R_{ij} A{j} S{j}
\\end{aligned}
```
with `C(PI)` is the observed counts in a detector channel `PI`, `T` is the observation time,
`ARF(E)` is the effective area of the instrument and `RMF(E, PI)` [we use `RMF(PI, E)`] is the unitless response matrix.
(Handbook of X-ray Astronomy Chapter 2011, Chapter 5 pg 87, by K. Arnaud, R. Smith and A. Siemiginowska)

This function takes the combined RMF and ARF as the response function. This is to recalculating it on every call.
Some people format the RMF as RMF(E, PI), others as RMF(PI, E). This latter convention is used by CIAO and by us.
"""
function apply_response_function(counts_per_bin::Vector{T}, response::Matrix{T}, exposure_time::T)::Vector{T} where {T<:AbstractFloat}
    @argcheck length(counts_per_bin) == size(response)[2]

    time_scaled_counts = counts_per_bin * exposure_time

    # All we have to do is matrix multiplication
    mult = response * time_scaled_counts
    return mult
end


function hydrogen_number_density(gas_density)
    # Abundances of elements, relative to nH. Taken from Anders & Grevesse (1989) https://doi.org/10.1016/0016-7037(89)90286-X
    abundance = 10 .^ ([12.00, 10.99, 8.56, 8.05, 8.93, 8.09, 6.33, 7.58, 6.47, 7.55, 7.21, 6.56, 6.36, 7.67, 6.25] .- 12)
    nucleon_total = [1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0, 28.0, 32.0, 40.0, 40.0, 56.0, 59.0]
    return gas_density / (m_p * dot(abundance, nucleon_total))
end

"""
    make_observation(temperature, density, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre, centre_radius, mask=nothing)

Generate an image of the cluster given functions for the radial profile of gas temperature and electron density and assorted observational parameters.

`shape` describes the size of the source number as a number of spatial bins in each dimension.
The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
This area is assumed to be square with the edge angle giving the side length.
The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
The response function includes both the RMF and ARF, as described in `apply_response_function`.
"""
function make_observation(
    temperature::Function,
    density::Function,
    z,
    shape,
    pixel_edge_angle::A,
    emission_model,
    exposure_time::T,
    response_function,
    centre::NTuple{2,<:DimensionfulAngles.Angle},
    centre_radius;
    mask::Union{Matrix{Bool},Nothing}=nothing,
    limit::Unitful.Length=Quantity(Inf, u"Mpc")
)::Array{Union{Float64,Missing},3} where {A<:DimensionfulAngles.Angle,T<:Unitful.Time}
    pixel_edge_length = ustrip(u"radᵃ", pixel_edge_angle) * angular_diameter_dist(cosmo, z)
    centre_length = ustrip.(u"radᵃ", centre) .* angular_diameter_dist(cosmo, z)
    radii_x, radii_y = shape ./ 2

    function radius_at_index(i, j, radii_x, radii_y, pixel_edge_length, centre_length)
        x = (i - radii_x) * pixel_edge_length - centre_length[1]
        y = (j - radii_y) * pixel_edge_length - centre_length[2]
        abs(hypot(x, y))
    end

    # Setting the min radius proportional to R500 throws off the results in favour of large mass
    # by varying the number of pixels contributing to the likelihood.
    min_radius = centre_radius * pixel_edge_length

    shortest_radius = min(radii_x * pixel_edge_length, radii_y * pixel_edge_length)
    if shortest_radius <= min_radius
        error("Minimum radius $min_radius greater than oberved radius in at least one direction ($shortest_radius).")
    end

    @mpirankeddebug "Creating brightness interpolation"
    brightness_radii = min_radius:pixel_edge_length:(hypot(radii_x + 1, radii_y + 1)*pixel_edge_length+hypot(centre_length...))

    brightness_line = [ustrip.(Float64, u"cm^(-2)/s", x) for x in surface_brightness.(
        brightness_radii,
        temperature,
        density,
        z,
        limit, # testing found that using 10Mpc as bound did not affect results
        Ref(emission_model),
        pixel_edge_angle
    )]

    brightness_interpolation = linear_interpolation(brightness_radii, brightness_line, extrapolation_bc=Throw())

    @mpirankeddebug "Calculating counts"
    resp = ustrip.(u"cm^2", response_function)
    exp_time = ustrip(u"s", exposure_time)
    counts = Array{Union{Float64,Missing}}(undef, size(response_function, 1), shape...)

    if isnothing(mask)
        mask = zeros(Bool, shape...)
    end

    for j in 1:shape[2]
        for i in 1:shape[1]
            radius = radius_at_index(i, j, radii_x, radii_y, pixel_edge_length, centre_length)
            if radius < min_radius
                counts[:, i, j] .= missing
            elseif mask[i, j]
                counts[:, i, j] .= missing
            else
                brightness = brightness_interpolation(radius)
                counts[:, i, j] .= apply_response_function(brightness, resp, exp_time)
            end
        end
    end

    return counts
end
"""
    make_observation(temperature, density, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre, centre_radius, mask=nothing)

Unitless wrapper for [make_observation](@ref)
"""
function make_observation(
    temperature::Function,
    density::Function,
    z,
    shape,
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