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
    pixel_edge_length
)
    @argcheck limit > 0u"Mpc"

    function integrand(l, params)
        s, temp = params
        r = hypot(s, l)
        kbT = ustrip(u"keV", temp(r))
        ρ = ustrip(u"cm^-3", density(r) / μ_e)

        f = model(kbT, ρ)

        @assert all(isfinite, f) "f with l=$l, s=$s (∴ r=$s, kbT=$kbT and ρ=$ρ) is $f"

        return f
    end

    problem = IntegralProblem(integrand, 0.0u"Mpc", limit, [projected_radius, temperature])
    sol = solve(problem, QuadGKJL(); reltol=1e-3, abstol=1e-3u"Mpc")

    @assert all(isfinite, sol.u)

    # display(sol.u)

    (1 / (4π * (1 + z)^4)) * (π^2 / (60^2 * 180^2)) * 2 * sol.u

end