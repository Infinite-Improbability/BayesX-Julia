export Model_Einasto

"""
    Model_Einasto(MT_200, fg_200, a, α, β, γ, c_500_GNFW, z)

Create functions for gas temperature and gas mass density as a function of radius.

Based on the Einasto mass density profile and the GNFW gas pressure profile.
"""
function Model_Einasto(
    MT_200::Unitful.Mass,
    fg_200,
    α,
    a,
    b,
    c,
    c_500_GNFW;
    z
)::NTuple{2,Function}
    # Move some parameters into a struct?

    @mpirankeddebug "Model" MT_200 fg_200

    @argcheck MT_200 > 0u"Msun"
    @argcheck 1 > fg_200 > 0
    @argcheck a > 0
    @argcheck c_500_GNFW > 0
    @argcheck (b - c_500_GNFW) > 0

    # Calculate NFW concentration parameter
    # This is equation 4 from Neto et al. 2007.
    # It assumes a relaxed halo and has different values in their full sample
    # Kinda sketch, I rather fit r200 or c200 as a prior
    c_200 = 5.26 * (MT_200 * cosmo.h / (10^14)u"Msun")^(-0.1) / (1 + z)
    # c_200_DM = 5.26 * (((MT_200 * cosmo.h) / 1e14)^(-0.1)) * (1 / (1 + z(k)))
    # Why does it have the redshift dependence?

    # Calculate gas mass
    Mg_200 = MT_200 * fg_200

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # And get R200 and NFW scale radius
    r_200 = uconvert(u"Mpc", cbrt((3 * MT_200) / (4π * 200 * ρ_crit_z)))
    rs = uconvert(u"Mpc", r_200 / c_200)

    # Calculate radius where mean density enclosed is@argcheck limit >
    # radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = ρ_crit_z * (200 / 3) * c_200^3 / (log(1 + c_200) - c_200 / (1 + c_200))

    r_500 = r_200 / 1.5

    # Set GNFW scale radius
    rp = uconvert(u"Mpc", r_500 / c_500_GNFW)

    function lower_gamma(s, x)
        try
            return gamma(s) - gamma(s, x)
        catch e
            @error("Gamma error, s cannot be over 51", s, x)
            throw(e)
        end
    end

    function pei_integrand(r, params)
        rs, α, rp, a, b, c = params
        r^3 * (b * (r / rp)^a + c) /
        lower_gamma(3 / α, 2 / α * (r / rs)^α) /
        (r / rp)^c / (1 + (r / rp)^a)^((a + b - c) / a)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    @mpirankeddebug "Integrating to find Pei"
    integral = IntegralProblem(pei_integrand, 0.0u"Mpc", r_200, (rs, α, rp, a, b, c))
    vol_int_200 = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4").u
    Pei_GNFW::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρ_s * rs^3 / α *
                                          (α / 2)^(3 / α) * exp(2 / α) *
                                          Mg_200 / vol_int_200
    @assert Pei_GNFW > 0u"Pa"
    @mpirankeddebug "Pei calculation complete"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{T})::Unitful.Density{T} where {T<:AbstractFloat}
        let
            ρ_s = ρ_s
            rs = rs
            rp = rp
            α = α
            a = a
            b = b
            c = c
            (μ_e / μ) * (1 / (4π * G)) * (Pei_GNFW / ρ_s) / rs^3 *
            α / (α / 2)^(3 / α) / exp(2 / α) *
            r / lower_gamma(3 / α, 2 / α * (r / rs)^α) *
            (r / rp)^(-c) * (1 + (r / rp)^a)^(-(a + b - c) / a) * (b * (r / rp)^a + c)
        end
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{T})::Unitful.Energy{T} where {T<:AbstractFloat}
        let
            ρ_s = ρ_s
            rs = rs
            rp = rp
            α = α
            a = a
            b = b
            c = c
            4π * μ * G * ρ_s * (rs^3) *
            (α / 2)^(3 / α) * exp(2 / α) / α *
            lower_gamma(3 / α, 2 / α * (r / rs)^α) / r *
            (1 + (r / rp)^α) / (b * (r / rp)^a + c)
        end
    end

    return gas_temperature, gas_density
end
"""
Mass is in solar masses.
"""
function Model_Einasto(
    MT_200::Real,
    fg_200,
    α,
    a,
    b,
    c,
    c_500_GNFW;
    z
)::NTuple{2,Function}
    Model_Einasto(
        MT_200 * 1u"Msun",
        fg_200,
        α,
        a,
        b,
        c,
        c_500_GNFW;
        z=z
    )
end