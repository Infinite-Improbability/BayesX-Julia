export Model_NFW

"""
    Model_NFW(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z, Δ=500)

Generate a cluster profile based on the NFW mass density and GNFW gas density profiles.

Uses the model from [olamaieSimpleParametricModel2012](@cite),
which is based on the NFW dark matter density profile and the GNFW gas pressure profile.

Keyword arguments are redshift and overdensity.

Returns functions for gas temperature and gas mass density as a function of radius.
"""
function Model_NFW(
    MT_Δ::Unitful.Mass,
    fg_Δ,
    c_Δ_dm,
    α,
    β,
    γ,
    c_Δ_GNFW;
    z,
    Δ=500
)::NTuple{2,Function}
    # Move some parameters into a struct?

    @mpirankeddebug "NFW" MT_Δ fg_Δ α β γ c_Δ_GNFW z
    @mpirankeddebug "NFW" MT_Δ fg_Δ α β γ c_Δ_GNFW z

    # Note the +1 in many likelihoods
    # This is so that something like fg_Δ=0 doesn't return a likelihood of zero despite failing the fg_Δ>0 constraint
    priorcheck(MT_Δ > 0u"Msun", -1e100(1 - ustrip(u"Msun", MT_Δ))) # MT_Δ is negative so we subtract it
    priorcheck(1 > fg_Δ > 0, -1e100(1 + abs(fg_Δ)))
    priorcheck(α > 0, -1e100(1 - α))
    priorcheck(c_Δ_GNFW > 0, -1e100(1 - c_Δ_GNFW))
    priorcheck(β > c_Δ_GNFW, -1e100(1 + (c_Δ_GNFW - β)))
    priorcheck(c_Δ_dm > 0, -1e100(1 - c_Δ_dm))

    # Calculate gas mass
    Mg_Δ = MT_Δ * fg_Δ

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get RΔ and NFW scale radius
    r_Δ = uconvert(u"Mpc", cbrt((3 * MT_Δ) / (4π * Δ * ρ_crit_z)))
    r_s = uconvert(u"Mpc", r_Δ / c_Δ_dm)

    # Calculate NFW characteristic overdensity
    ρ_s = ρ_crit_z * (Δ / 3) * c_Δ_dm^3 / (log(1 + c_Δ_dm) - c_Δ_dm / (1 + c_Δ_dm))

    # Set GNFW scale radius
    r_p = uconvert(u"Mpc", r_Δ / c_Δ_GNFW)

    # Some helper functions
    """The radius dependent part of the gas density function"""
    function gnfw_gas_radial_term(
        r::Unitful.Length{<:Real},
        r_s::Unitful.Length{<:Real}, # NFW
        r_p::Unitful.Length{<:Real}, # GNFW
        α,
        β,
        γ
    )::Unitful.Length{Float64}
        r = uconvert(u"Mpc", r)
        r_s = uconvert(u"Mpc", r_s)
        r_p = uconvert(u"Mpc", r_p)
        r / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
        (r / r_p)^(-γ) *
        (1 + (r / r_p)^α)^(-(α + β - γ) / α) *
        (β * (r / r_p)^α + γ)
    end

    """An integral over radius that is equal to the gas
    density to a proportionality constant"""
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{<:Real},
        r_s::Unitful.Length{<:Real}, # NFW
        r_p::Unitful.Length{<:Real}, # GNFW
        α,
        β,
        γ
    )::Unitful.Volume{Float64}
        s = r^2 * gnfw_gas_radial_term(r, r_s, r_p, α, β, γ)

        @assert isfinite(s) "Not finite with $r, $r_s, $r_p, $α, $β, $γ"

        return s
    end
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{<:Real},
        p::Tuple
    )
        gnfw_gas_mass_integrand(r, p...)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    integral = IntegralProblem(
        gnfw_gas_mass_integrand,
        0.0u"Mpc",
        r_Δ,
        (r_s, r_p, α, β, γ)
    )
    vol_int_Δ = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4").u
    Pei_GNFW::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_Δ / vol_int_Δ
    @assert Pei_GNFW > 0u"Pa"
    @mpirankeddebug "Pei calculation complete"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{<:Real})::Unitful.Density{Float64}
        r = abs(r)
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
    function gas_temperature(r::Unitful.Length{<:Real})::Unitful.Energy{Float64}
        r = abs(r)
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
    Model_NFW(MT_Δ::Real, fg_Δ, α, β, γ, c_Δ_GNFW; z, Δ=500)

Unitless wrapper for [`Model_NFW`](@ref). Mass is in solar masses.
"""
function Model_NFW(
    MT_Δ::Real,
    fg_Δ,
    c_Δ_dm,
    α,
    β,
    γ,
    c_Δ_GNFW;
    z,
    Δ=500
)::NTuple{2,Function}
    Model_NFW(
        MT_Δ * 1u"Msun",
        fg_Δ,
        c_Δ_dm,
        α,
        β,
        γ,
        c_Δ_GNFW;
        z=z,
        Δ=Δ
    )
end