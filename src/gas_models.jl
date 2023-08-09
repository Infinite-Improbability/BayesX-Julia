using ArgCheck
using Unitful, UnitfulAstro, DimensionfulAngles
using PhysicalConstants.CODATA2018: G
using Integrals
# using Optimization
# using OptimizationOptimJL

include("params.jl")
include("emission.jl")


"""
    p_crit(z)

Calculate the critical density at some redshift `z`.
"""
ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)


"""
    Model_NFW_GNFW(MT_200, fg_200, α, β, γ, c_500_GNFW, z)

Create functions for gas temperature and electron density as a function of radius.

Uses the model from Olamaie 2012 (doi:10.1111/j.1365-2966.2012.20980.x),
which is based on the NFW dark matter density profile and the GNFW gas pressure profile.
"""
function Model_NFW_GNFW(
    MT_200::Unitful.Mass,
    fg_200,
    α,
    β,
    γ,
    c_500_GNFW,
    z,
)::NTuple{2,Function}
    # Move some parameters into a struct?

    @mpirankeddebug "Model" MT_200 fg_200

    @argcheck MT_200 > 0u"Msun"
    @argcheck 1 > fg_200 > 0
    @argcheck α > 0
    @argcheck c_500_GNFW > 0
    @argcheck (β - c_500_GNFW) > 0

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

    return gas_temperature, gas_density
end
"""
Mass is in solar masses.
"""
function Model_NFW_GNFW(
    MT_200::Real,
    fg_200,
    α,
    β,
    γ,
    c_500_GNFW,
    z,
)::NTuple{2,Function}
    Model_NFW_GNFW(
        MT_200 * 1u"Msun",
        fg_200,
        α,
        β,
        γ,
        c_500_GNFW,
        z
    )
end