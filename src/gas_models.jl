using ArgCheck
using Unitful, UnitfulAstro, UnitfulAngles
using PhysicalConstants.CODATA2018: G
using Integrals


include("params.jl")
include("emission.jl")
include("mpi.jl")


"""
    p_crit(z)

Calculate the critical density at some redshift `z`.
"""
ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)


"""
    Model_NFW_GNFW_GNFW(MT_200, fg_200, a_GNFW, b_GNFW, c_GNFW, c_500_GNFW, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function)

Calculate predicted counts using a physical model based NFW-GNFW profiles as described in Olamaie 2012.

The first six parameters are model priors. 
`z` is redshift and `shape` describes the size of the source number as a number of spatial bins in each dimension.
The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
This area is assumed to be square with the edge angle giving the side length.
The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
The response function includes both the RMF and ARF, as described in `apply_response_function`.
"""
function Model_NFW_GNFW(
    MT_200::T,
    fg_200::T,
    a_GNFW::T,
    b_GNFW::T,
    c_GNFW::T,
    c_500_GNFW::T,
    z::T,
    shape::Vector{N},
    pixel_edge_angle::Quantity{T,NoDims},
    emission_model,
    exposure_time::Unitful.Time,
    response_function::Matrix,
)::Array{Float64} where {N<:Integer,T<:AbstractFloat}
    # Move some parameters into an object?

    @mpirankeddebug "Model called with parameters MT_200=$MT_200, fg_200=$fg_200"

    @argcheck MT_200 > 0
    @argcheck fg_200 > 0
    @argcheck a_GNFW > 0
    # @argcheck c_500_GNFW > 0
    # @argcheck (b_GNFW - c_500_GNFW) > 0

    MT_200 *= 1u"Msun" # todo: unitful as primary method, with wrapper to add units?

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

    # Sketchy way to get R500 from r200
    # TODO: Do better
    r_500 = r_200 / 1.5
    # c_500 = r_500 / r_s

    # Set GNFW scale radius
    r_p = uconvert(u"Mpc", r_500 / c_500_GNFW)

    # Some helper functions
    """The radius dependent part of the gas density function"""
    function gnfw_gas_radial_term(
        r::Unitful.Length{T},
        r_s::Unitful.Length{T}, # NFW
        r_p::Unitful.Length{T}, # GNFW
        a,
        b,
        c
    ) where {T<:AbstractFloat}
        r / (log(1 + r / r_s) - (1 + r_s / r)^(-1)) *
        (r / r_p)^(-c) *
        (1 + (r / r_p)^a)^(-(a + b - c) / a) *
        (b * (r / r_p)^a + c)
    end

    """An integral over radius that is equal to the gas
    density to a proportionality constant"""
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{T},
        r_s::Unitful.Length{T}, # NFW
        r_p::Unitful.Length{T}, # GNFW
        a,
        b,
        c
    ) where {T<:AbstractFloat}
        s = r^2 * gnfw_gas_radial_term(r, r_s, r_p, a, b, c)

        @assert isfinite(s) "Not finite with $r, $r_s, $r_p, $a, $b, $c"

        return s
    end
    function gnfw_gas_mass_integrand(
        r::Unitful.Length{T},
        p::Vector{Quantity{T}}
    ) where {T<:AbstractFloat}
        gnfw_gas_mass_integrand(r, p...)
    end

    # Calculate Pei, normalisation coefficent for GNFW pressure

    @mpirankeddebug "Integrating to find Pei"
    integral = IntegralProblem(
        gnfw_gas_mass_integrand,
        0.0u"Mpc",
        r_200,
        [r_s, r_p, a_GNFW, b_GNFW, c_GNFW]
    )
    vol_int_200 = solve(integral, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc^4")
    Pei_GNFW = (μ / μ_e) * G * ρ_s * r_s^3 * Mg_200_DM / vol_int_200.u
    @mpirankeddebug "Pei calculation complete"

    @assert Pei_GNFW > 0u"Pa"

    """Calculate gas density at some radius"""
    function gas_density(r::Unitful.Length)::Unitful.Density
        (μ_e / μ) * (1 / (4π * G)) *
        (Pei_GNFW / ρ_s) * (1 / r_s^3) *
        gnfw_gas_radial_term(r, r_s, r_p, a_GNFW, b_GNFW, c_GNFW)
    end

    """Calculate gas temperature at some radius"""
    function gas_temperature(r::Unitful.Length)::Unitful.Energy
        4π * μ * G * ρ_s * (r_s^3) *
        ((log(1 + r / r_s) - (1 + r_s / r)^(-1)) / r) *
        (1 + (r / r_p)^a_GNFW) * (b_GNFW * (r / r_p)^a_GNFW + c_GNFW)^(-1)
    end

    # Calculate source brightness at various points
    # TODO: Moving center
    pixel_edge_length = ustrip(u"rad", pixel_edge_angle) * angular_diameter_dist(cosmo, z)
    radii_x, radii_y = ceil.(Int64, shape ./ 2)

    radius_at_coords(x, y) = hypot(x, y) * pixel_edge_length

    radius_at_cell = Matrix{typeof(0.0u"Mpc")}(undef, radii_x, radii_y)

    # TODO: Are we transposing?
    for y in 1:radii_y
        for x in 1:radii_x
            radius_at_cell[x, y] = radius_at_coords(x, y)
        end
    end

    @mpirankeddebug "Calculating brightness"

    brightness = surface_brightness.(
        radius_at_cell,
        gas_temperature,
        gas_density, # TODO: Calculate nH instead of using ne
        z,
        20 * max(radii_x, radii_y) * pixel_edge_length,
        Ref(emission_model),
        pixel_edge_length
    )

    counts = Matrix{Vector{Float64}}(undef, size(brightness)...)

    @mpirankeddebug "Applying response function"
    # Stripping units here is an attempt to fix performance issues
    # I think the root cause is using generic instead of specialised matrix multiplication
    # We may be able to clean that up with better typing instead of stripping unitd
    resp = ustrip.(u"cm^2", response_function)
    exp_time = ustrip(u"s", exposure_time)
    for i in eachindex(brightness)
        counts[i] = apply_response_function(ustrip.(u"cm^-2/s", brightness[i]), resp, exp_time)
    end
    #

    # Potential optimisations
    # Supply integrals as Vector
    # Eliminate duplicate radii

    # TODO: Conversion from rate to counts
    # For now we handwave it with a fixed factor

    # Debug code for verifing matrix completion and radius calculation
    # for y in 1:radii_y
    #     for x in 1:radii_x
    #         counts[x, y] = repeat([ustrip(u"Mpc", radius_at_coords(x, y))], 26)
    #     end
    # end

    return complete_matrix(counts, shape)
end
function Model_NFW_GNFW(
    MT_200::Unitful.Mass,
    fg_200,
    a_GNFW,
    b_GNFW,
    c_GNFW,
    c_500_GNFW,
    z,
    shape::Vector{N},
    pixel_edge_angle::Quantity{T,NoDims},
    emission_model,
    exposure_time::Unitful.Time,
    response_function::Matrix,
)::Array{Float64} where {N<:Integer,T<:AbstractFloat}
    Model_NFW_GNFW(
        ustrip(u"Msun", MT_200),
        fg_200,
        a_GNFW,
        b_GNFW,
        c_GNFW,
        c_500_GNFW,
        z,
        shape,
        pixel_edge_angle,
        emission_model,
        exposure_time,
        response_function
    )
end

"""
    complete_matrix(m::Matrix, shape::Vector)

Complete matrix expands predictions from 1/4 of the sky to the whole sky.

Takes a small matrix `m` of vectors of counts per channel and treats it as the quadrants of a larger matrix of size `shape`.
This is expanded by channels to create a 3D array of counts for `(channel, x, y)`.
"""
function complete_matrix(m::Matrix, shape::Vector{N})::Array{Float64} where {N<:Int}
    @mpirankeddebug "Completing matrix"
    new = Array{Float64}(undef, length(m[1]), shape...)

    # increasing row is increasing x
    # increasing column is increasing y

    radii = ceil.(Int64, shape / 2)

    # top left
    for y in 1:radii[2]
        for x in 1:radii[1]
            new[:, radii[1]+1-x, radii[2]+1-y] = m[x, y]
        end
    end

    # bottom left
    for y in 1:radii[2]
        for x in 1:radii[1]
            new[:, radii[1]+x, radii[2]+1-y] = m[x, y]
        end
    end

    # top right
    for y in 1:radii[2]
        for x in 1:radii[1]
            new[:, radii[1]+1-x, radii[2]+y] = m[x, y]
        end
    end

    # bottom right
    for y in 1:radii[2]
        for x in 1:radii[1]
            new[:, radii[1]+x, radii[2]+y] = m[x, y]
        end
    end

    @assert all(isfinite, new)

    @mpirankeddebug "Matrix reshaped"
    return new
end
