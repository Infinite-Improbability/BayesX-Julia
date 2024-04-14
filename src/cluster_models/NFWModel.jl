export NFWModel


"""
    NFWModel(ρs, rs, rp, α, β, γ, P0)

Generate a cluster profile based on the NFW mass density and GNFW gas density profiles.

Uses the model from [olamaieSimpleParametricModel2012](@cite),
which is based on the NFW dark matter density profile and the GNFW gas pressure profile.
"""
struct NFWModel{L<:Unitful.Length} <: ClusterModel
    ρs
    rs::L
    rp::L
    α
    β
    γ
    P0
end
# Normalise units on init

"""The radius dependent part of the gas density function"""
function gnfw_radial_term(r, rs, rp, α, β, γ)
    r /
    (log1p(r / rs) - (1 + rs / r)^(-1)) *
    (r / rp)^(-γ) *
    (1 + (r / rp)^α)^(-(α + β - γ) / α) *
    (β * (r / rp)^α + γ)
end
gnfw_radial_term(cm::NFWModel, r::Unitful.Length{<:Real})::Unitful.Length{Float64} = gnfw_radial_term(r, cm.rs, cm.rp, cm.α, cm.β, cm.γ)

"""
function gnfw_mass_integrand(cm::NFWModel, r::Unitful.Length{<:Real})::Unitful.Volume{Float64}

An integral over radius that is equal to the gas
density to a proportionality constant"""
function gnfw_mass_integrand(r, rs, rp, α, β, γ)::Unitful.Volume{Float64}

    s = r^2 * gnfw_radial_term(r, rs, rp, α, β, γ)
    @assert isfinite(s) "s is not finite with r = $r, rs = $rs, rp = $rp, α = $α, β = $β, γ = $γ"

    return s
end
gnfw_mass_integrand(r, p) = gnfw_mass_integrand(r, p...)

"""
    NFWModel(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z, Δ=500)

Constructor for [`NFWModel`](@ref) using standard cluster properties.

`MT_Δ` is assumed to be in solar masses if not specified. Keyword arguments are redshift and overdensity.
"""
function NFWModel(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z, Δ=500)
    # @mpirankeddebug "NFW" MT_Δ fg_Δ α β γ c_Δ_GNFW z Δ

    # Note the +1 in many likelihoods
    # This is so that something like fg_Δ=0 doesn't return a likelihood of zero despite failing the fg_Δ>0 constraint
    priorcheck(MT_Δ > 0u"Msun", -1e100(1 - ustrip(u"Msun", MT_Δ))) # MT_Δ is negative so we subtract it
    priorcheck(1 > fg_Δ > 0, -1e100(1 + abs(fg_Δ)))
    priorcheck(α > 0, -1e100(1 - α))
    priorcheck(c_Δ_GNFW > 0, -1e100(1 - c_Δ_GNFW))
    priorcheck(β > γ, -1e100(1 + (γ - β)))
    priorcheck(c_Δ_dm > 0, -1e100(1 - c_Δ_dm))

    # Calculate gas mass at overdensity Δ
    Mg_Δ = MT_Δ * fg_Δ

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # Get overdensity and scale radii
    rΔ = uconvert(u"kpc", cbrt((3 * MT_Δ) / (4π * Δ * ρ_crit_z)))
    rs = uconvert(u"kpc", rΔ / c_Δ_dm) # NFW scale radius
    rp = uconvert(u"kpc", rΔ / c_Δ_GNFW) # GNFW scale radius

    # Calculate NFW characteristic overdensity
    ρs = ρ_crit_z * (Δ / 3) * c_Δ_dm^3 / (log1p(c_Δ_dm) - c_Δ_dm / (1 + c_Δ_dm))

    # Calculate P0, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find P0"
    ifunc = IntegralFunction(gnfw_mass_integrand)
    integral = IntegralProblem(
        ifunc,
        (0.0u"kpc", rΔ),
        (rs, rp, α, β, γ)
    )
    vol_int_Δ = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1u"kpc^4").u
    P0::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρs * rs^3 * Mg_Δ / vol_int_Δ
    @assert isfinite(P0)
    @assert P0 > 0u"Msun/kpc/s^2"
    @mpirankeddebug "P0 calculation complete" uconvert(u"Msun/kpc/s^2", P0) ρ_s rs Mg_Δ vol_int_Δ

    # TODO: Iterative corrections to estimate that include gas density (see Javid et al. 2019)

    return NFWModel(ρs, rs, rp, α, β, γ, P0)
end
NFWModel(MT_Δ::Real, fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z, Δ=500) = NFWModel(MT_Δ * 1u"Msun", fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z=z, Δ=Δ)

function density(cm::NFWModel, r::Unitful.Length{<:Real})
    r = abs(r)

    (μ_e / μ) * (1 / (4π * G)) *
    (cm.P0 / cm.ρs) * (1 / cm.rs^3) *
    gnfw_radial_term(cm, r)
end

function temperature(cm::NFWModel, r::Unitful.Length{<:Real})
    r = abs(r)

    4π * μ * G * cm.ρs * (cm.rs^3) *
    ((log1p(r / cm.rs) - (1 + cm.rs / r)^(-1)) / r) *
    (1 + (r / cm.rp)^cm.α) *
    (cm.β * (r / cm.rp)^cm.α + cm.γ)^(-1)
end
