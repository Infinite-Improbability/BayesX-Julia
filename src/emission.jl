using Integrals

include("mekal.jl")

"""Observed surface brightness"""
function surface_brightness(
    projected_radius::Unitful.Length,
    temperature::Function,
    density::Function,
    z,
    limit::Unitful.Length,
    model,
    pixel_edge_length::Unitful.Length,
    exposure_time::Unitful.Time
)
    @argcheck limit > 0u"Mpc"

    function integrand(l, params)
        s, temp = params
        r = hypot(s, l)

        # should be passing n_h not n_e to integral
        # and using n_e n_h instead of n_e^2
        f = model(ustrip(u"keV", temp(r)), ustrip(u"cm^-3", density(r) / μ_e)) * 1u"cm^3/s" * (density(r) / μ_e)^2

        # TODO: Switch to in place

        @assert all(isfinite, f) "f with l=$l, s=$s (∴ r=$s, kbT=$kbT and ρ=$ρ) is $f"

        return f
    end

    problem = IntegralProblem(integrand, 0.0u"Mpc", limit, [projected_radius, temperature])
    sol = solve(problem, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc/cm^2/s")

    @assert all(isfinite, sol.u)

    # doubling solution to account for integral bounds
    # applying exposure area and time
    # applying XSPEC normalisation
    # the n_e n_h density normalisation has been applied in the integrand, approximating n_e n_h as n_e^2
    # this is not ideal, should use the explicit value
    # I don't know where XSPEC gets 10^-14 so I'm just trusting their docs
    2 * sol.u * pixel_edge_length^2 * exposure_time * 10^-14 / 4π / (1 + z)^2 / angular_diameter_dist(cosmo, z)^2

end