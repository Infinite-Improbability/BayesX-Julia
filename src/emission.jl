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
    exposure_time::Unitful.Time=1u"s"
)::Vector{Float64}
    @argcheck limit > 0u"Mpc"

    function integrand(l, params)
        s, temp = params
        r = hypot(s, l)

        f = model(ustrip(u"keV", temp(r)), ustrip(u"cm^-3", density(r) / μ_e)) * 1u"cm^(-2)/s"

        # TODO: Switch to in place

        @assert all(isfinite, f) "f with l=$l, s=$s (∴ r=$s, kbT=$kbT and ρ=$ρ) is $f"

        return f
    end

    problem = IntegralProblem(integrand, 0.0u"Mpc", limit, [projected_radius, temperature])
    sol = solve(problem, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc/cm^2/s")

    @assert all(isfinite, sol.u)

    # display(sol.u)

    # x2 to account for half integral bounds
    (1 / (4π * (1 + z)^4)) * (π^2 / (60^2 * 180^2)) * 2 * sol.u * pixel_edge_length^2 * exposure_time

end