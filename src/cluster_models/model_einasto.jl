using SpecialFunctions

export Model_Einasto

function lower_gamma(s, x)
    try
        return gamma(s) - gamma(s, x)
    catch e
        @error("Gamma error, s cannot be over 51", s, x)
        throw(e)
    end
end

function safe_lower_gamma(s, x)
    g = lower_gamma(s, x)
    return g < 0 ? 0.0 : g
end

"""
    Model_Einasto(MT_200, fg_200, α, a, b, c, c_500_GNFW, z)

Generate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.

Returns functions for gas temperature and gas mass density as a function of radius.
See Appendix A, Olamaie et al. 2015 (doi:10.1093/mnras/stu2146)

If `α>2` then the gas temperature starts increasing to physically improbable
levels (orders exceeding 10¹¹ keV). We thus constrain the value to below that. 
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

    @mpirankeddebug "Einasto" MT_200 fg_200 α a b c c_500_GNFW z


    # Note the +1 in many likelihoods
    # This is so that something like fg_200=0 doesn't return a likelihood of zero
    priorcheck(MT_200 > 0u"Msun", -1e100(1 - ustrip(u"Msun", MT_200))) # MT_200 is negative so we subtract it
    priorcheck(1 > fg_200 > 0, -1e100(1 + abs(fg_200)))
    priorcheck(a > 0, -1e100(1 - a))
    priorcheck(c_500_GNFW > 0, -1e100(1 - c_500_GNFW))
    priorcheck(b > c_500_GNFW, -1e100(1 + (c_500_GNFW - b)))
    priorcheck(α < 2, -1e100α) # α > 2 leads to temperature going up past some r
    priorcheck(α >= 0.06, -1e100(1 - α)) # if 3/α > 51 then lower_gamma(3/α, x) throws an error

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

    function pei_integrand(r, params)
        rs, α, rp, a, b, c = params
        r^3 * (b * (r / rp)^a + c) /
        safe_lower_gamma(3 / α, 2 / α * (r / rs)^α) /
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

    @assert ρ_s > 0u"kg/m^3"
    @assert rs > 0u"Mpc"
    @assert rp > 0u"Mpc"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{T})::Unitful.Density{T} where {T<:AbstractFloat}
        # lower_gamma(3/1.5, 0.5e-31 < 0 presumably occurs elsewhere
        # this is a problem
        r = abs(r)
        res = let
            ρ_s = ρ_s
            rs = rs
            rp = rp
            α = α
            a = a
            b = b
            c = c
            (μ_e / μ) * (1 / (4π * G)) * (Pei_GNFW / ρ_s) / rs^3 *
            α / (α / 2)^(3 / α) / exp(2 / α) *
            r / safe_lower_gamma(3 / α, 2 / α * (r / rs)^α) *
            (r / rp)^(-c) * (1 + (r / rp)^a)^(-(a + b - c) / a) * (b * (r / rp)^a + c)
        end
        isfinite(res) ? res : 0.0u"kg/cm^3"
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{T})::Unitful.Energy{T} where {T<:AbstractFloat}
        r = abs(r)
        res = let
            ρ_s = ρ_s
            rs = rs
            rp = rp
            α = α
            a = a
            b = b
            c = c
            4π * μ * G * ρ_s * (rs^3) *
            (α / 2)^(3 / α) * exp(2 / α) / α *
            safe_lower_gamma(3 / α, 2 / α * (r / rs)^α) / r *
            (1 + (r / rp)^α) / (b * (r / rp)^a + c)
        end
        isfinite(res) ? res : 0.0u"keV"
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