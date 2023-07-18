using Integrals
using LinearAlgebra: dot

include("mekal.jl")

"""
    surface_brightness(projected_radius, temperature, density, z, limit, model, pixel_edge_length)

Calculate the observed surface_brightness at some projected radius on the sky.

The temperature and density functions should take the true distance from origin (cluster center) and return the appropriate value, with units.
The limit controls limits on LOS radius during integration. Theoretically it should be ±infinity, but we may approximate it as a finite value.
The pixel edge length is the distance observed on the sky, determined from the pixel_edge_angle at source redshift.
"""
function surface_brightness(
    projected_radius::Unitful.Length,
    temperature::Function,
    density::Function,
    z,
    limit::Unitful.Length,
    model,
    pixel_edge_length::Unitful.Length,
)
    @argcheck limit > 0u"Mpc"

    function integrand(l, params)
        s, temp = params
        r = hypot(s, l)

        nH = hydrogen_number_density(density(r))

        f = model(ustrip(u"keV", temp(r)), ustrip(u"cm^-3", nH)) * 1u"cm^3/s" * (density(r) / μ_e) * nH

        # TODO: Switch to in place

        @assert all(isfinite, f) "f with l=$l, s=$s (∴ r=$s, kbT=$kbT and ρ=$ρ) is $f"

        return f
    end

    problem = IntegralProblem(integrand, 0.0u"Mpc", limit, [projected_radius, temperature])
    sol = solve(problem, QuadGKJL(); reltol=1e-3, abstol=1e-3u"cm^-2/s")

    @assert all(isfinite, sol.u)

    # doubling solution to account for integral bounds
    # applying exposure area and time
    # applying XSPEC normalisation
    # the n_e n_h density normalisation has been applied in the integrand, approximating n_e n_h as n_e^2
    # this is not ideal, should use the explicit value
    # I don't know where XSPEC gets 10^-14 so I'm just trusting their docs
    2 * sol.u * pixel_edge_length^2 * 10^-14 / 4π / (1 + z)^2 / angular_diameter_dist(cosmo, z)^2
end

"""
    apply_response_function(counts_per_bin::Vector, response::Matrix, exposure_time::Unitful.Time)::Vector{Float64}

Applies the response function ``RSP(PI, E) = RMF(PI, E) ∘ ARF(E)`` to all energy bins and returns adjusted counts per bin.
```math
\\begin{aligned}
C(PI) &= T \\int RMF(PI, E) ⋅ ARF(E)sum ⋅ S(E) dE \\
&≈ T \\sum_{j} R_{ij} A{j} S{j}
\\end{aligned}
```
with `C(PI)`` is the observed counts in a detector channel `PI`, `T` is the observation time,
`ARF(E)` is the effective area of the instrument and `RMF(E, PI)` is the unitless response matrix.
(Handbook of X-ray Astronomy Chapter 2011, Chapter 5 pg 87, by K. Arnaud, R. Smith and A. Siemiginowska)

This function takes the combined RMF and ARF as the response function. This is to recalculating it on every call.
Some people format the RMF as RMF(PI, E). This convention is used by CIAO, for example.
"""
function apply_response_function(counts_per_bin::Vector, response::Matrix, exposure_time)::Vector{Float64}
    # Argcheck would really be better here but we want it to skip it in high performance situations
    # display(counts_per_bin)
    # display(response)
    @assert length(counts_per_bin) == size(response)[2] "There are $(length(counts_per_bin)) energy bins but the response matrix has $(size(response, 2)) columns"

    # display(response)
    # display(counts_per_bin)
    # display(exposure_time)


    time_scaled_counts = counts_per_bin * exposure_time

    # All we have to do is matrix multiplication
    mult = response * time_scaled_counts
    return mult
end

function hydrogen_number_density(gas_density)
    abundance = 10 .^ ([12.00, 10.99, 8.56, 8.05, 8.93, 8.09, 6.33, 7.58, 6.47, 7.55, 7.21, 6.56, 6.36, 7.67, 6.25] .- 12)
    nucleon_total = [1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0, 28.0, 32.0, 40.0, 40.0, 56.0, 59.0]
    return gas_density / (m_p * dot(abundance, nucleon_total))
end