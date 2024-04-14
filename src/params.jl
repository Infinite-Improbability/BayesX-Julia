using Cosmology
using Unitful, DimensionfulAngles, UnitfulAstro
using PhysicalConstants.CODATA2018: m_p, G, k_B

"""
Cosmology assumed by BayesX
"""
const cosmo = cosmology(
    h=0.6774,
    Neff=3.046,
    OmegaK=0,
    OmegaM=0.3089,
    OmegaR=nothing,
    Tcmb=2.7255,
    w0=-1,
    wa=0) # todo: make configurable

"""Mean mass per gas particle = 0.6mₚ"""
const μ = 0.6 * m_p

"""Mean gas mass per electron = 1.14mₚ"""
const μ_e = 1.14 * m_p

function angle_to_length(θ::A, redshift) where {A<:DimensionfulAngles.Angle}
    return ustrip(u"radᵃ", θ) * angular_diameter_dist(cosmo, redshift)
end

function get_centre_indices(centre_x::A, centre_y::A, pixel_edge_angle::A, shape::NTuple{2,<:Integer}) where {A<:DimensionfulAngles.Angle}
    radii_x, radii_y = shape ./ 2

    i = centre_x / pixel_edge_angle + radii_x
    j = centre_y / pixel_edge_angle + radii_y

    return (i, j)
end

@derived_dimension SurfaceDensity Unitful.𝐋^-2
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension NumberDensityRate Unitful.𝐋^-3 / Unitful.𝐓