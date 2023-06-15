using ArgCheck
using Unitful, UnitfulAstro, UnitfulAngles
using PhysicalConstants.CODATA2018: G
using Integrals
using SpectralFitting
using ProgressMeter
using Interpolations


include("params.jl")
include("utils.jl")


ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)

"""Observed surface brightness"""
function surface_brightness(
    projected_radius::Unitful.Length,
    temperature::Function,
    z,
    limit::Unitful.Length,
    model
)
    @argcheck limit > 0u"Mpc"

    function integrand(l, p)
        s, temp = p
        r::Unitful.Length = hypot(s, l)
        kbT = ustrip(u"keV", temp(r))

        # TODO: Better emission model

        f = model(kbT)

        @assert !any(isnan, f)

        return f
    end

    # TODO: Try infinite bounds
    problem = IntegralProblem(integrand, 0.0u"Mpc", limit, [projected_radius, temperature])
    sol = solve(problem, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc")

    (1 / (4π * (1 + z)^4)) * (π^2 / (60^2 * 180^2)) * 2 * ustrip.(u"Mpc", sol.u)

end


"""Calculate predicted counts using a physical model based NFW-GNFW profiles
as described in Olamaie 2012.
"""
function Model_NFW_GNFW(;
    MT_200::Unitful.Mass,
    fg_200,
    a_GNFW,
    b_GNFW,
    c_GNFW,
    c_500_GNFW,
    z,
    radius_steps::Integer=10,
    # radius_limits::Tuple{Unitful.Length},
    energy_limits::Vector{T},
    n_energy_bins::Integer,
    shape::Vector{N}=[12, 12],
    pixel_edge_angle
) where {T<:Unitful.Energy,N<:Integer}
    # Move some parameters into an object?

    @argcheck MT_200 > 0u"Msun"
    @argcheck fg_200 > 0
    @argcheck a_GNFW > 0
    # @argcheck c_500_GNFW > 0
    # @argcheck (b_GNFW - c_500_GNFW) > 0

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

    # And get R200 and NFW scale radius
    r_200 = uconvert(u"Mpc", cbrt((3 * MT_200) / (4π * 200 * ρ_crit_z)))
    r_s = uconvert(u"Mpc", r_200 / c_200)

    # Calculate radius where mean density enclosed is@argcheck limit >
    # radii = LogRange(radius_limits..., radius_steps)

    # Calculate NFW characteristic overdensity
    ρ_s = ρ_crit_z * (200 / 3) * c_200 / (log(1 + c_200) - c_200 / (1 + c_200))

    # Sketchy way to get R500 from r200
    # TODO: Do better
    r_500 = r_200 / 1.5
    # c_500 = r_500 / r_s

    # Set GNFW scale radius
    r_p = uconvert(u"Mpc", r_500 / c_500_GNFW)

    # Some helper functions

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
        r^2 * gnfw_gas_radial_term(r, r_s, r_p, a, b, c)
    end
    function gnfw_gas_mass_integrand(
        r::Unitful.Length,
        p
    )
        gnfw_gas_mass_integrand(r, p...)
    end


    """The radius dependent part of the gas density function"""
    function gnfw_gas_radial_term(
        r::Unitful.Length,
        r_s::Unitful.Length, # NFW
        r_p::Unitful.Length, # GNFW
        a,
        b,
        c
    )
        r / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
        (r / r_p)^(-c) *
        (1 + (r / r_p)^a)^(-(a + b - c) / a) *
        (b * (r / r_p)^a + c)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure
    # Find source or verify this equation

    @info "Integrating to find Pei"
    integral = IntegralProblem(
        gnfw_gas_mass_integrand,
        0.0u"Mpc",
        r_200,
        [r_s, r_p, a_GNFW, b_GNFW, c_GNFW]
    )
    vol_int_200 = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4")
    Pei_GNFW = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_200_DM / vol_int_200.u
    @info "Pei calculation complete"

    @assert Pei_GNFW > 0u"Pa"

    # Coeffient to integral for surface brightness in Olamaie 2015 Eq8
    # Called xfluxsec1 in BayesX-Fortran
    # sx_coefficent = 1 / (4π * (1 + z)^4)

    # Create energy bins
    @argcheck energy_limits[1] < energy_limits[2]
    @argcheck energy_limits[1] >= 0u"keV"

    energy_bins = LinRange(energy_limits..., n_energy_bins)

    # Calculate absorption
    # TODO: Confirm this is transmission
    # absorption_model = PhotoelectricAbsorption()
    # transmission = invokemodel(energy_bins, absorption_model)

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length)::Unitful.Density
        (μ_e / μ) * (1 / (4π * G)) *
        (Pei_GNFW / ρ_s) * (1 / r_s^3) *
        gnfw_gas_radial_term(r_200, r_s, r_p, a_GNFW, b_GNFW, c_GNFW)
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length)::Unitful.Energy
        4π * μ * G * ρ_s * (r_s^3) *
        ((log(1 + r / r_s) - (1 + r_s / r)^(-1)) / r) *
        (1 + (r / r_p)^a_GNFW) * (b_GNFW * (r / r_p)^a_GNFW + c_GNFW)^(-1)
    end

    # ρg_200::Unitful.Density = gas_density(r_200)
    # Tg_200::Unitful.Temperature = gas_temperature(r_200)

    # Calculate source brightness at various points
    # TODO: Moving center
    pixel_edge_length = ustrip(u"rad", pixel_edge_angle) * angular_diameter_dist(cosmo, z)
    radii_x, radii_y = ceil.(Int64, shape ./ 2)

    radius_at_cell = zeros(Float64, (radii_x, radii_y))
    counts = zeros(Float64, (radii_x, radii_y))

    # TODO: Are we transposing?
    for y in 1:radii_y
        for x in 1:radii_x
            radius_at_cell[x, y] = hypot(x, y)
        end
    end
    radius_at_cell *= pixel_edge_length

    @info "Preparing model"
    model(kbT) = PhotoelectricAbsorption() * (XS_BremsStrahlung(T=FitParam(kbT)) + BlackBody(kT=FitParam(kbT)))
    temps = 0:0.01:3
    surrogate = linear_interpolation(temps, invokemodel.(Ref(ustrip.(u"keV", energy_bins)), model.(temps)))

    @info "Generating counts"

    counts = surface_brightness.(
        radius_at_cell,
        gas_temperature,
        z,
        20 * max(radii_x, radii_y) * pixel_edge_length,
        Ref(surrogate)
    )
    @info "Done"

    # Potential optimisations
    # Supply integrals as Vector
    # Eliminate duplicate radii

    # return counts
end

# function xray_flux_coefficent()

# end

# @xspecmodel :C_mekal struct XS_Mekal{T,F} <: AbstractSpectralModel{T,Additive}
#     "Normalisation"
#     K::T
#     "Plasma Temperature"
#     t::T
#     "Hydrogen Density"
#     ρ::T
#     "Metal Abundances"
#     a::T
#     "Redshift"
#     z::T
#     "Switch"
#     s::Int
# end
# function XS_Mekal(;
#     K=FitParam(1.0),
#     t=FitParam(8.0),
#     ρ=FitParam(10.0),
#     a=FitParam(1.0),
#     z=FitParam(0.1),
#     s=0
# )
#     XS_Mekal{typeof(K),SpectralFitting.FreeParameters{(:K, :t)}}(
#         K, t, ρ, a, z, s
#     )
# end
# SpectralFitting.register_model_data(XS_Mekal, "mekal1.dat")

# model = XS_Mekal(t=FitParam(8.0), ρ=FitParam(12.0), z=FitParam(0.1))
# invokemodel(collect(0.1:0.1:2), model)

@time Model_NFW_GNFW(
    MT_200=5e14u"Msun",
    fg_200=0.13,
    a_GNFW=1.0620,
    b_GNFW=5.4807,
    c_GNFW=0.3292,
    c_500_GNFW=1.156,
    z=0.1,
    radius_steps=10,
    energy_limits=[0.3u"keV", 3u"keV"],
    n_energy_bins=27,
    shape=[12, 12],
    pixel_edge_angle=0.492u"arcsecond"
)

@time Model_NFW_GNFW(
    MT_200=5e14u"Msun",
    fg_200=0.13,
    a_GNFW=1.0620,
    b_GNFW=5.4807,
    c_GNFW=0.3292,
    c_500_GNFW=1.156,
    z=0.1,
    radius_steps=10,
    energy_limits=[0.3u"keV", 3u"keV"],
    n_energy_bins=27,
    shape=[12, 12],
    pixel_edge_angle=0.492u"arcsecond"
)