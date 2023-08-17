export Model_NFW

"""
    Model_NFW(MT_200::Unitful.Mass, fg_200, a, b, c, c_500_GNFW, z)

Generate a cluster profile based on the NFW mass density and GNFW gas density profiles.

Uses the model from Olamaie 2012 (doi:10.1111/j.1365-2966.2012.20980.x),
which is based on the NFW dark matter density profile and the GNFW gas pressure profile.
Returns functions for gas temperature and gas mass density as a function of radius.
"""
function Model_NFW(
    MT_200::Unitful.Mass,
    fg_200,
    a,
    b,
    c,
    c_500_GNFW;
    z
)::NTuple{2,Function}
    # Move some parameters into a struct?

    @mpirankeddebug "NFW" MT_200 fg_200 a b c c_500_GNFW z

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
        r::Unitful.Length{<:Real},
        r_s::Unitful.Length{<:Real}, # NFW
        r_p::Unitful.Length{<:Real}, # GNFW
        a,
        b,
        c
    )::Unitful.Length{Float64}
        r / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
        (r / r_p)^(-c) *
        (1 + (r / r_p)^a)^(-(a + b - c) / a) *
        (b * (r / r_p)^a + c)
    end

    """An integral over radius that is equal to the gas
    density to a proportionality constant"""
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{<:Real},
        r_s::Unitful.Length{<:Real}, # NFW
        r_p::Unitful.Length{<:Real}, # GNFW
        a,
        b,
        c
    )::Unitful.Volume{Float64}
        s = r^2 * gnfw_gas_radial_term(r, r_s, r_p, a, b, c)

        @assert isfinite(s) "Not finite with $r, $r_s, $r_p, $a, $b, $c"

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
        r_200,
        (r_s, r_p, a, b, c)
    )
    vol_int_200 = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4").u
    Pei_GNFW::Unitful.Pressure{Float64} = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_200_DM / vol_int_200
    @assert Pei_GNFW > 0u"Pa"
    @mpirankeddebug "Pei calculation complete"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length{<:Real})::Unitful.Density{Float64}
        r = abs(r)
        let
            ρ_s = ρ_s
            r_s = r_s
            r_p = r_p
            a = a
            b = b
            c = c
            (μ_e / μ) * (1 / (4π * G)) *
            (Pei_GNFW / ρ_s) * (1 / r_s^3) *
            gnfw_gas_radial_term(r, r_s, r_p, a, b, c)
        end
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length{<:Real})::Unitful.Energy{Float64}
        r = abs(r)
        let
            ρ_s = ρ_s
            r_s = r_s
            r_p = r_p
            a = a
            b = b
            c = c
            4π * μ * G * ρ_s * (r_s^3) *
            ((log(1 + r / r_s) - (1 + r_s / r)^(-1)) / r) *
            (1 + (r / r_p)^a) * (b * (r / r_p)^a + c)^(-1)
        end
    end

    return gas_temperature, gas_density
end
"""
    Model_NFW(MT_200::Real, fg_200, a, b, c, c_500_GNFW, z)

Generate a cluster profile based on the NFW mass density and GNFW gas density profiles.

Mass is in solar masses.
"""
function Model_NFW(
    MT_200::Real,
    fg_200,
    a,
    b,
    c,
    c_500_GNFW;
    z
)::NTuple{2,Function}
    Model_NFW(
        MT_200 * 1u"Msun",
        fg_200,
        a,
        b,
        c,
        c_500_GNFW;
        z=z
    )
end