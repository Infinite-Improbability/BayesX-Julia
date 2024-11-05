using SpecialFunctions

export Model_Einasto

# SpecialFunctions has several algorithms avaliable for this - there may be a better choice
lower_gamma(s, x) = gamma(s) * gamma_inc(s, x)[1]

"""
    Model_Einasto(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z, Δ=500)

Generate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.

See Appendix A of [olamaieBAYESXBayesianInference2015](@cite) for a derivation
of the model based on these profiles.

Returns functions for gas temperature and gas mass density as a function of radius.

If `1/n>2` then the gas temperature starts increasing to physically improbable
levels (orders exceeding 10¹¹ keV).
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
    # priorcheck(n  > 0.5, -1e100/n) # 1/n > 2 leads to temperature going up past some r
    priorcheck(3n <= 51, -1e100(3n - 15)) # if 3n > 51 then lower_gamma(3n, x) throws an error
    priorcheck(c_Δ_dm > 0, -1e100(1 - c_Δ_dm))


    # Calculate gas mass
    Mg_Δ = MT_Δ * fg_Δ

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get RΔ and NFW scale radius
    r_Δ = uconvert(u"kpc", cbrt((3 * MT_Δ) / (4π * Δ * ρ_crit_z)))
    r_s = uconvert(u"kpc", r_Δ / c_Δ_dm)

    # Calculate radius where mean density enclosed is@argcheck limit >
    # radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = uconvert(u"Msun/kpc^3", ρ_crit_z * (Δ / 3) * c_Δ_dm^3 / (log1p(c_Δ_dm) - c_Δ_dm / (1 + c_Δ_dm)))

    # Set GNFW scale radius
    r_p = uconvert(u"kpc", r_Δ / c_Δ_GNFW)

    function pei_integrand(r, params)
        rs, n, rp, α, β, γ = params
        r^3 * (β * (r / rp)^α + γ) /
        lower_gamma(3n, 2n * (r / rs)^(1 / n)) /
        (r / rp)^γ / (1 + (r / rp)^α)^((α + β - γ) / α)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    ifunc = IntegralFunction(pei_integrand)
    integral = IntegralProblem(ifunc, (0.0u"kpc", r_Δ), (r_s, n, r_p, α, β, γ))
    vol_int_Δ = try
        solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"kpc^4").u
    catch e
        if e isa DomainError
            @mpirankeddebug "Pei integral failed, returning -1e-100"
            throw(PriorError(-1e-100))
        else
            rethrow()
        end
    end

    Pei_GNFW = (μ / μ_e) * G * ρ_s * r_s^3 * n * (1 / 2n)^(3n) * exp(2n) * Mg_Δ / vol_int_Δ
    @assert isfinite(Pei_GNFW)
    @assert Pei_GNFW > 0u"Msun/kpc/s^2"

    @assert ρ_s > 0u"Msun/kpc^3"
    @assert r_s > 0u"kpc"
    @assert r_p > 0u"kpc"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{T})::Unitful.Density{T} where {T<:AbstractFloat}
        r = abs(r)
        res = let ρ_s = ρ_s, rs = r_s, rp = r_p, n = n, α = α, β = β, γ = γ
            (μ_e / μ) * (1 / (4π * G)) * (Pei_GNFW / ρ_s) / rs^3 /
            n / (1 / 2n)^(3n) / exp(2n) *
            r / lower_gamma(3n, 2n * ustrip(r / rs)^(1 / n)) *
            (r / rp)^(-γ) * (1 + (r / rp)^α)^(-(α + β - γ) / α) * (β * (r / rp)^α + γ)
        end
        @assert isfinite(res)
        return res
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{T})::Unitful.Energy{T} where {T<:AbstractFloat}
        r = abs(r)
        res = let ρ_s = ρ_s, rs = r_s, rp = r_p, n = n, α = α, β = β, γ = γ
            4π * μ * G * ρ_s * (rs^3) *
            (1 / 2n)^(3n) * exp(2n) * n *
            lower_gamma(3n, 2n * (r / rs)^(1 / n)) / r *
            (1 + (r / rp)^(1 / n)) / (β * (r / rp)^α + γ)
        end
        @assert isfinite(res)
        return res
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