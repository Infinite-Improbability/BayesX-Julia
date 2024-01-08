using SpecialFunctions

export Model_Einasto

function lower_gamma(s, x)
    try
        return gamma(s) - gamma(s, x)
    catch e
        @error("Gamma error, s cannot be over 51", s, x)
        rethrow()
    end
end

function safe_lower_gamma(s, x)
    g = lower_gamma(s, x)
    return g < 0 ? 0.0 : g
end

"""
    Model_Einasto(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z, Δ=500)

Generate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.

See Appendix A of [olamaieBAYESXBayesianInference2015](@cite) for a derivation
of the model based on these profiles.

Returns functions for gas temperature and gas mass density as a function of radius.

If `n>2` then the gas temperature starts increasing to physically improbable
levels (orders exceeding 10¹¹ keV). We thus constrain the value to below that. 
"""
function Model_Einasto(
    MT_Δ::Unitful.Mass,
    fg_Δ,
    c_Δ_dm,
    n,
    α,
    β,
    γ,
    c_Δ_GNFW;
    z,
    Δ=500
)::NTuple{2,Function}
    # Move some parameters into a struct?

    @mpirankeddebug "Einasto" MT_Δ fg_Δ n α β γ c_Δ_GNFW z

    # Note the +1 in many likelihoods
    # This is so that something like fg_Δ=0 doesn't return a likelihood of zero
    priorcheck(MT_Δ > 0u"Msun", -1e100(1 - ustrip(u"Msun", MT_Δ))) # MT_Δ is negative so we subtract it
    priorcheck(1 > fg_Δ > 0, -1e100(1 + abs(fg_Δ)))
    priorcheck(α > 0, -1e100(1 - α))
    priorcheck(c_Δ_GNFW > 0, -1e100(1 - c_Δ_GNFW))
    priorcheck(β > c_Δ_GNFW, -1e100(1 + (c_Δ_GNFW - β)))
    priorcheck(n < 2, -1e100n) # n > 2 leads to temperature going up past some r
    priorcheck(n >= 0.06, -1e100(1 - n)) # if 3/n > 51 then lower_gamma(3/n, x) throws an error
    priorcheck(c_Δ_dm > 0, -1e100(1 - c_Δ_dm))


    # Calculate gas mass
    Mg_Δ = MT_Δ * fg_Δ

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get RΔ and NFW scale radius
    r_Δ = uconvert(u"Mpc", cbrt((3 * MT_Δ) / (4π * Δ * ρ_crit_z)))
    r_s = uconvert(u"Mpc", r_Δ / c_Δ_dm)

    # Calculate radius where mean density enclosed is@argcheck limit >
    # radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = ρ_crit_z * (Δ / 3) * c_Δ_dm^3 / (log(1 + c_Δ_dm) - c_Δ_dm / (1 + c_Δ_dm))

    # Set GNFW scale radius
    r_p = uconvert(u"Mpc", r_Δ / c_Δ_GNFW)

    function pei_integrand(r, params)
        rs, n, rp, α, β, γ = params
        r^3 * (β * (r / rp)^α + γ) /
        safe_lower_gamma(3 / n, 2 / n * (r / rs)^n) /
        (r / rp)^γ / (1 + (r / rp)^α)^((α + β - γ) / α)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    integral = IntegralProblem(pei_integrand, 0.0u"Mpc", r_Δ, (r_s, n, r_p, α, β, γ))
    vol_int_Δ = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4").u
    Pei_GNFW::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρ_s * r_s^3 / n *
                                          (n / 2)^(3 / n) * exp(2 / n) *
                                          Mg_Δ / vol_int_Δ
    @assert Pei_GNFW > 0u"Pa"
    @mpirankeddebug "Pei calculation complete"

    @assert ρ_s > 0u"kg/m^3"
    @assert r_s > 0u"Mpc"
    @assert r_p > 0u"Mpc"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{T})::Unitful.Density{T} where {T<:AbstractFloat}
        # lower_gamma(3/1.5, 0.5e-31) < 0
        # presumably occurs elsewhere
        # this is a problem
        r = abs(r)
        res = let ρ_s = ρ_s, rs = r_s, rp = r_p, n = n, α = α, β = β, γ = γ
            (μ_e / μ) * (1 / (4π * G)) * (Pei_GNFW / ρ_s) / rs^3 *
            n / (n / 2)^(3 / n) / exp(2 / n) *
            r / safe_lower_gamma(3 / n, 2 / n * (r / rs)^n) *
            (r / rp)^(-γ) * (1 + (r / rp)^α)^(-(α + β - γ) / α) * (β * (r / rp)^α + γ)
        end
        isfinite(res) ? res : 0.0u"kg/cm^3"
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{T})::Unitful.Energy{T} where {T<:AbstractFloat}
        r = abs(r)
        res = let ρ_s = ρ_s, rs = r_s, rp = r_p, n = n, α = α, β = β, γ = γ
            4π * μ * G * ρ_s * (rs^3) *
            (n / 2)^(3 / n) * exp(2 / n) / n *
            safe_lower_gamma(3 / n, 2 / n * (r / rs)^n) / r *
            (1 + (r / rp)^n) / (β * (r / rp)^α + γ)
        end
        isfinite(res) ? res : 0.0u"keV"
    end

    return gas_temperature, gas_density
end
"""
    Model_Einasto(MT_Δ::Real, fg_Δ, n, α, β, γ, c_Δ_GNFW; z, Δ=500)

Unitless wrapper for [`Model_Einasto`](@ref). Mass is in solar masses.
"""
function Model_Einasto(
    MT_Δ::Real,
    fg_Δ,
    c_Δ_dm,
    n,
    α,
    β,
    γ,
    c_Δ_GNFW;
    z,
    Δ=500
)::NTuple{2,Function}
    Model_Einasto(
        MT_Δ * 1u"Msun",
        fg_Δ,
        c_Δ_dm,
        n,
        α,
        β,
        γ,
        c_Δ_GNFW;
        z=z,
        Δ=Δ
    )
end