using ArgCheck
using Unitful, UnitfulAstro
using PhysicalConstants.CODATA2018: G
using Integrals
using SpectralFitting: PhotoelectricAbsorption

include("params.jl")
include("utils.jl")

ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)

"""An integral over radius that is equal to the gas
    density to a proportionality constant"""
function gnfw_gas_mass_integrand(
    r::Unitful.Length,
    r_s::Unitful.Length, # NFW
    r_p::Unitful.Length, # GNFW
    a,
    b,
    c
)
    r^3 / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
    (r / r_p)^(-c) *
    (1 + (r / r_p)^a)^(-(a + b - c) / a) *
    (b * (r / r_p)^a + c)
end
function gnfw_gas_mass_integrand(
    r::Unitful.Length,
    p
)
    gnfw_gas_mass_integrand(r, p...)
end

"""Calculate predicted counts using a physical model based NFW-GNFW profiles
as described in Olamaie 2012.
"""
function Model_NFW_GNFW(
    MT_200::Unitful.Mass,
    fg_200,
    a_GNFW,
    b_GNFW,
    c_GNFW,
    c_500_GNFW,
    z;
    radius_steps <: Integer=10,
    radius_limits::Tuple{Unitful.Length{T}},
    energy_limits::Tuple{Unitful.Energy},
    energy_bins <: Integer
)
    # Move some parameters into an object?

    @argcheck MT_200 > 0
    @argcheck fg_200 > 0
    @argcheck a_GNFW > 0
    @argcheck >= 0

    # Calculate NFW concentration parameter
    # This is equation 4 from Neto et al. 2007.
    # It assumes a relaxed halo and has different values in their full sample
    # Kinda sketch, I rather fit r200 or c200 as a prior
    c_200 = 5.26 * (MT_200 * cosmo.h / (10^14)u"Msun")^(-0.1)
    # c_200_DM = 5.26 * (((MT_200 * cosmo.h) / 1e14)^(-0.1)) * (1 / (1 + z(k)))
    # Why does BayesX have the redshift dependence?

    # Calculate gas mass
    Mg_200_DM = MT_200 * fg_200

    # Calculate critical density at current redshift
    ρ_crit_z = ρ_crit(z)

    # Calculate radius where mean density enclosed is 200 times the critical density
    r_200 = cbrt((3 * MT_200) / (4π * 200 * ρ_crit_z))

    # Calculate NFW scale radius
    r_s = r_200 / c_200

    # Discretize r
    @argcheck radius_limits[1] < radius_limits[2]
    @argcheck radius_limits[1] >= 0u"Mpc"
    radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = (200 / 3) * c_200 / (log(1 + c_200) - c_200 / (1 + c_200))

    # Sketchy way to get R500 from r200
    # TODO: Do better
    r_500 = r_200 / 1.5
    c_500 = r_500 / r_s

    # Set GNFW scale radius
    r_p = r_500 / c_500_GNFW

    # Calculate Pei, normalisation coefficent for GNFW pressure
    # Find source or verify this equation
    integral = IntegralProblem(gnfw_gas_mass_integrand, 0, r_200, p=(r_s, ρ_s, r_p, a, b, c))
    vol_int_200 = solve(integral, QuadGKJL).u
    Pei_GNFW::Unitful.Pressure = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_200_DM / vol_int_200

    @assert Pei_GNFW > 0

    # Coeffient to integral for surface brightness in Olamaie 2015 Eq8
    # Called xfluxsec1 in BayesX-Fortran
    sx_coefficent = 1 / (4π * (1 + z)^4)

    # Create energy bins
    @argcheck energy_limits[1] < energy_limits[2]
    @argcheck energy_limits[1] >= 0u"keV"

    energy_bins = LinRange(energy_limits..., energy_bins)

    # DO ABSORPTION CALCULATIONS
    model = PhotoelectricAbsorption()

    # DO XRAY_FLUX_COEFF CALCULATIONS

end

function xray_flux_coefficent()

end