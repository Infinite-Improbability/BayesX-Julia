using ArgCheck
using Unitful, UnitfulAstro, DimensionfulAngles
using PhysicalConstants.CODATA2018: G
using Integrals
using Optimization
using OptimizationOptimJL

include("params.jl")
include("emission.jl")

export Model_NFW_GNFW_Priors

"""
    p_crit(z)

Calculate the critical density at some redshift `z`.
"""
ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)


"""
    Model_NFW_GNFW_GNFW(; MT_200, fg_200, α, β, γ, c_500_GNFW, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre_x, centre_y, centre_radius)

Calculate predicted counts using a physical model based NFW-GNFW profiles as described in Olamaie 2012.

The first six parameters are model priors. 
`z` is redshift and `shape` describes the size of the source number as a number of spatial bins in each dimension.
The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
This area is assumed to be square with the edge angle giving the side length.
The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
The response function includes both the RMF and ARF, as described in `apply_response_function`.
"""
function Model_NFW_GNFW(;
    MT_200::Unitful.Mass{<:Real},
    fg_200,
    α,
    β,
    γ,
    c_500_GNFW,
    z,
    shape::AbstractVector{<:Integer},
    pixel_edge_angle::DimensionfulAngles.Angle{<:Real},
    emission_model,
    exposure_time::Unitful.Time{<:Real},
    response_function::Matrix,
    centre_x::DimensionfulAngles.Angle{<:Real},
    centre_y::DimensionfulAngles.Angle{<:Real},
    centre_radius
)::Array{Float64}
    # Move some parameters into a struct?

    @mpirankeddebug "Model" MT_200 fg_200

    @argcheck MT_200 > 0u"Msun"
    @argcheck 1 > fg_200 > 0
    @argcheck α > 0
    @argcheck c_500_GNFW > 0
    @argcheck (β - c_500_GNFW) > 0

    centre_coordinates = (centre_x, centre_y)

    # Calculate NFW concentration parameter
    # This is equation 4 from Neto et al. 2007.
    # It assumes a relaxed halo and has different values in their full sample
    # Kinda sketch, I rather fit r200 or c200 as a prior
    c_200 = 5.26 * (MT_200 * cosmo.h / (10^14)u"Msun")^(-0.1) / (1 + z)
    # c_200_DM = 5.26 * (((MT_200 * cosmo.h) / 1e14)^(-0.1)) * (1 / (1 + z(k)))
    # Why does it have the redshift dependence?

    # Calculate gas mass
    Mg_200_DM = MT_200 * fg_200

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get R200 and NFW scale radius
    r_200 = uconvert(u"Mpc", cbrt((3 * MT_200) / (4π * 200 * ρ_crit_z)))
    r_s = uconvert(u"Mpc", r_200 / c_200)

    # Calculate radius where mean density enclosed is@argcheck limit >
    # radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = ρ_crit_z * (200 / 3) * c_200^3 / (log(1 + c_200) - c_200 / (1 + c_200))

    # function m(r, p)
    #     x = r[1] * 1u"Mpc"
    #     return abs(ustrip(u"Msun", 4π * ρ_s * r_s^3 * (log(1 + x / r_s) - (1 + r_s / x)^(-1)) - 4π / 3 * x^3 * p[1] * ρ_crit_z))
    # end

    # opf = OptimizationFunction(m, AutoForwardDiff())
    # op = OptimizationProblem(opf, [ustrip(u"Mpc", r_200 / 1.5)], [500], lb=[0], ub=[ustrip(u"Mpc", r_200)])
    # r_500 = solve(op, GradientDescent()).u[1] * 1u"Mpc"
    r_500 = r_200 / 1.5

    # Set GNFW scale radius
    r_p = uconvert(u"Mpc", r_500 / c_500_GNFW)

    # Some helper functions
    """The radius dependent part of the gas density function"""
    function gnfw_gas_radial_term(
        r::Unitful.Length{T},
        r_s::Unitful.Length{T}, # NFW
        r_p::Unitful.Length{T}, # GNFW
        α,
        β,
        γ
    )::Unitful.Length{T} where {T<:AbstractFloat}
        r / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
        (r / r_p)^(-γ) *
        (1 + (r / r_p)^α)^(-(α + β - γ) / α) *
        (β * (r / r_p)^α + γ)
    end

    """An integral over radius that is equal to the gas
    density to a proportionality constant"""
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{T},
        r_s::Unitful.Length{T}, # NFW
        r_p::Unitful.Length{T}, # GNFW
        α,
        β,
        γ
    )::Unitful.Volume{T} where {T<:AbstractFloat}
        s = r^2 * gnfw_gas_radial_term(r, r_s, r_p, α, β, γ)

        @assert isfinite(s) "Not finite with $r, $r_s, $r_p, $α, $β, $γ"

        return s
    end
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{T},
        p::Vector{Quantity{T}}
    ) where {T<:AbstractFloat}
        gnfw_gas_mass_integrand(r, p...)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    integral = IntegralProblem(
        gnfw_gas_mass_integrand,
        0.0u"Mpc",
        r_200,
        [r_s, r_p, α, β, γ]
    )
    vol_int_200 = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4").u
    Pei_GNFW::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_200_DM / vol_int_200
    @assert Pei_GNFW > 0u"Pa"
    @mpirankeddebug "Pei calculation complete"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{T})::Unitful.Density{T} where {T<:AbstractFloat}
        let
            ρ_s = ρ_s
            r_s = r_s
            r_p = r_p
            α = α
            β = β
            γ = γ
            (μ_e / μ) * (1 / (4π * G)) *
            (Pei_GNFW / ρ_s) * (1 / r_s^3) *
            gnfw_gas_radial_term(r, r_s, r_p, α, β, γ)
        end
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{T})::Unitful.Energy{T} where {T<:AbstractFloat}
        let
            ρ_s = ρ_s
            r_s = r_s
            r_p = r_p
            α = α
            β = β
            γ = γ
            4π * μ * G * ρ_s * (r_s^3) *
            ((log(1 + r / r_s) - (1 + r_s / r)^(-1)) / r) *
            (1 + (r / r_p)^α) * (β * (r / r_p)^α + γ)^(-1)
        end
    end

    # Calculate source brightness at various points
    # TODO: Moving center
    pixel_edge_length = ustrip(u"radᵃ", pixel_edge_angle) * angular_diameter_dist(cosmo, z)
    centre_length = ustrip.(u"radᵃ", centre_coordinates) .* angular_diameter_dist(cosmo, z)
    radii_x, radii_y = shape ./ 2

    function radius_at_index(i, j, radii_x, radii_y, pixel_edge_length, centre_length)
        x = (i - radii_x) * pixel_edge_length - centre_length[1]
        y = (j - radii_y) * pixel_edge_length - centre_length[2]
        abs(hypot(x, y))
    end

    # min_radius = r_500 * 0.1
    min_radius = centre_radius * pixel_edge_length

    shortest_radius = min(radii_x * pixel_edge_length, radii_y * pixel_edge_length)
    if shortest_radius <= min_radius
        error("Minimum radius $min_radius greater than oberved radius in at least one direction ($shortest_radius).")
    end

    @mpirankeddebug "Creating brightness interpolation"
    brightness_radii = min_radius:pixel_edge_length:(hypot(radii_x + 1, radii_y + 1)*pixel_edge_length+hypot(centre_length...))
    brightness_line = [ustrip.(Float64, u"cm^(-2)/s", x) for x in surface_brightness.(
        brightness_radii,
        gas_temperature,
        gas_density,
        z,
        Quantity(Inf, u"Mpc"),
        Ref(emission_model),
        pixel_edge_angle
    )]
    brightness_interpolation = linear_interpolation(brightness_radii, brightness_line, extrapolation_bc=Throw())

    @mpirankeddebug "Calculating counts"
    resp = ustrip.(u"cm^2", response_function)
    exp_time = ustrip(u"s", exposure_time)
    counts = Array{Float64}(undef, size(resp, 1), shape...)

    for j in 1:shape[2]
        for i in 1:shape[1]
            radius = radius_at_index(i, j, radii_x, radii_y, pixel_edge_length, centre_length)
            if radius < min_radius
                counts[:, i, j] .= NaN
            else
                # if hydrogen_number_density(gas_density(radius)) > 1u"cm^-3"
                #     @error "Extreme gas density at" radius min_radius r_500 pixel_edge_angle pixel_edge_length
                # end
                brightness = brightness_interpolation(radius)
                counts[:, i, j] .= apply_response_function(brightness, resp, exp_time)
            end
        end
    end

    # replace!(i -> i < 0 ? 0 : i, counts)

    # @assert all(i -> i >= 0, counts)

    return counts
end
# function Model_NFW_GNFW(;
#     MT_200,
#     centre_x,
#     centre_y,
#     kwargs...
# )::Array{Float64}
#     display(kwargs)
#     Model_NFW_GNFW(
#         MT_200=MT_200 * 1u"Msun", # values passed in as priors may lack units
#         centre_x=centre_x .* 1u"arcsecondᵃ",
#         centre_y=centre_y .* 1u"arcsecondᵃ",
#         kwargs...
#     )
# end

struct Model_NFW_GNFW_Priors{M<:Prior,F<:Prior,A<:Prior,B<:Prior,G<:Prior,C<:Prior,X<:Prior,Y<:Prior} <: ModelPriors
    MT_200::M
    fg_200::F
    α::A
    β::B
    γ::G
    c_500_GNFW::C
    x::X
    y::Y
end