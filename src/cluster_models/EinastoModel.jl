using SpecialFunctions

export EinastoModel

# SpecialFunctions has several algorithms avaliable for this - there may be a better choice
lower_gamma(s, x) = gamma(s) * gamma_inc(s, x)[1]

"""
    EinastoModel(ρs, n, rs, rp, α, β, γ, P0)

Generate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.

See Appendix A of [olamaieBAYESXBayesianInference2015](@cite) for a derivation
of the model based on these profiles.
"""
struct EinastoModel{L<:Unitful.Length} <: ClusterModel
    ρs
    n
    rs::L
    rp::L
    α
    β
    γ
    P0
end
# Normalise units on init

function pei_integrand(r, rs, n, rp, α, β, γ)
    r^3 * (β * (r / rp)^α + γ) /
    lower_gamma(3n, 2n * ustrip(u"kpc/kpc", r / rs)^(1 / n)) /
    (r / rp)^γ / (1 + (r / rp)^α)^((α + β - γ) / α)
end
function pei_integrand(r, params)
    pei_integrand(r, params...)
end

"""
    EinastoModel(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z, Δ=500)

Constructor for [`EinastoModel`](@ref) using standard cluster properties.

`MT_Δ` is assumed to be in solar masses if not specified. Keyword arguments are redshift and overdensity.
"""
function EinastoModel(
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
)
    # Move some parameters into a struct?

    @mpirankeddebug "Einasto" MT_Δ fg_Δ n α β γ c_Δ_GNFW z

    # Note the +1 in many likelihoods
    # This is so that something like fg_Δ=0 doesn't return a likelihood of zero
    priorcheck(MT_Δ > 0u"Msun", -1e100(1 - ustrip(u"Msun", MT_Δ))) # MT_Δ is negative so we subtract it
    priorcheck(1 > fg_Δ > 0, -1e100(1 + abs(fg_Δ)))
    priorcheck(α > 0, -1e100(1 - α))
    priorcheck(c_Δ_GNFW > 0, -1e100(1 - c_Δ_GNFW))
    priorcheck(β > c_Δ_GNFW, -1e100(1 + (c_Δ_GNFW - β)))
    # priorcheck(n  > 0.5, -1e100/n) # 1/n > 2 leads to temperature going up to improbable levels past some r
    priorcheck(3n <= 51, -1e100(3n - 15)) # if 3n > 51 then lower_gamma(3n, x) throws an error
    priorcheck(c_Δ_dm > 0, -1e100(1 - c_Δ_dm))

    # Calculate gas mass
    Mg_Δ = MT_Δ * fg_Δ

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get RΔ and NFW scale radius
    rΔ = uconvert(u"kpc", cbrt((3 * MT_Δ) / (4π * Δ * ρ_crit_z)))
    rs = uconvert(u"kpc", rΔ / c_Δ_dm)

    # Calculate NFW characteristic overdensity
    ρs = uconvert(u"Msun/kpc^3", ρ_crit_z * (Δ / 3) * c_Δ_dm^3 / (log1p(c_Δ_dm) - c_Δ_dm / (1 + c_Δ_dm)))

    # Set GNFW scale radius
    rp = uconvert(u"kpc", rΔ / c_Δ_GNFW)

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    ifunc = IntegralFunction(pei_integrand)
    integral = IntegralProblem(ifunc, (0.0u"kpc", rΔ), (rs, n, rp, α, β, γ))
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

    P0 = (μ / μ_e) * G * ρs * rs^3 * n * (1 / 2n)^(3n) * exp(2n) * Mg_Δ / vol_int_Δ
    @assert isfinite(P0)
    @assert P0 > 0u"Msun/kpc/s^2"
    @mpirankeddebug "Pei calculation complete"

    @assert ρs > 0u"Msun/kpc^3"
    @assert rs > 0u"kpc"
    @assert rp > 0u"kpc"

    return EinastoModel(ρs, n, rs, rp, α, β, γ, uconvert(u"Msun/kpc/s^2", P0))
end
EinastoModel(MT_Δ::Real, fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z, Δ=500) = EinastoModel(MT_Δ * 1u"Msun", fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z=z, Δ=Δ)

gamma_term(cm::EinastoModel, r) = lower_gamma(3cm.n, 2cm.n * ustrip(u"kpc/kpc", r / cm.rs)^(1 / cm.n))


function density(cm::EinastoModel, r::Unitful.Length{<:Real})
    r = abs(r)

    res = (μ_e / μ) * (1 / (4π * G)) * (cm.P0 / cm.ρs) / cm.n / (1 / 2cm.n)^(3cm.n) / cm.rs^3 / exp(2cm.n) *
          r / gamma_term(cm, r) *
          (r / cm.rp)^(-cm.γ) * (1 + (r / cm.rp)^cm.α)^(-(cm.α + cm.β - cm.γ) / cm.α) * (cm.β * (r / cm.rp)^cm.α + cm.γ)
    @assert isfinite(res)

    return res
end

function temperature(cm::EinastoModel, r::Unitful.Length{<:Real})
    r = abs(r)

    res = 4π * μ * G * cm.ρs * (cm.rs^3) *
          (1 / 2cm.n)^(3cm.n) * exp(2cm.n) * cm.n *
          gamma_term(cm, r) / r *
          (1 + (r / cm.rp)^(1 / cm.n)) /
          (cm.β * (r / cm.rp)^cm.α + cm.γ)

    @assert isfinite(res)

    return res
end