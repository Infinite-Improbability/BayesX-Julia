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

    projected_radius = ustrip(u"Mpc", projected_radius)
    pixel_edge_length = ustrip(u"Mpc", pixel_edge_length)
    limit = ustrip(u"Mpc", limit)

    function integrand(vars, params)
        s, temp = params
        l, x, y = vars
        r = hypot(x + y, l) * 1u"Mpc"
        kbT = ustrip(u"keV", temp(r))
        ρ = ustrip(u"cm^-3", density(r) / μ_e)

        # TODO: Better emission model

        f = model(kbT, ρ)

        @assert all(isfinite, f) "f with l=$l, s=$s (∴ r=$s, kbT=$kbT and ρ=$ρ) is $f"

        return f
    end

    # TODO: Try infinite bounds
    upper_bounds = [Inf, hypot(projected_radius, pixel_edge_length), hypot(projected_radius, pixel_edge_length)]
    lower_bounds = [-Inf, projected_radius, projected_radius]

    problem = IntegralProblem(integrand, lower_bounds, upper_bounds, [projected_radius, temperature])
    sol = solve(problem, HCubatureJL(); reltol=1e-3, abstol=1e-3)

    @assert all(isfinite, sol.u)

    # display(sol.u)

    (1 / (4π * (1 + z)^4)) * (π^2 / (60^2 * 180^2)) * 2 * sol.u

end