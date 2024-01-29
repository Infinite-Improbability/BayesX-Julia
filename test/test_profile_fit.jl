using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using MPI
using LinearAlgebra: I

z = 0.1
shape = (5, 5)

pixel_edge_angle = 20.0u"arcsecondᵃ"
energy_bins = range(0.7u"keV", 7.0u"keV", step=0.01u"keV")
exposure_time = 3.0e6u"s"
response_function = 250u"cm^2" * Matrix(I, (length(energy_bins) - 1, length(energy_bins) - 1))
centre_radius = 0
integration_limit = 10.0u"Mpc"

emission_model = BayesJ.prepare_model_mekal(
    2.0e20u"cm^-2",
    energy_bins,
    z,
)

temperature, density = Model_NFW(3.e14, 0.13, 3.0, 1.0510, 5.4905, 0.3081, 1.177, z=z)

predicted_count_rate = BayesJ.make_observation(
    temperature,
    density,
    z,
    shape,
    pixel_edge_angle,
    emission_model,
    exposure_time,
    response_function,
    (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
    centre_radius,
    limit=integration_limit,
)

@assert all(isfinite, predicted_count_rate)

observation = pois_rand.(predicted_count_rate)
background_rate = rand()
background = Array{Int64}(undef, size(observation)...)
for i in eachindex(background)
    background[i] = pois_rand(background_rate) + 1
end

priors = [
    DeltaPrior("x0", 0.0),
    DeltaPrior("y0", 0.0),
    DeltaPrior("r0", 0.0), LogUniformPrior("ρ0", 1.e-28, 1.e-23), UniformPrior("T0", 0.0, 5.0),
    DeltaPrior("r1", 300.0), LogUniformPrior("ρ1", 1.e-28, 1.e-23), UniformPrior("T1", 0.0, 5.0),
    DeltaPrior("r2", 500.0), LogUniformPrior("ρ2", 1.e-28, 1.e-23), UniformPrior("T2", 0.0, 5.0),
    DeltaPrior("r3", 1000.0), LogUniformPrior("ρ3", 1.e-28, 1.e-23), UniformPrior("T3", 0.0, 5.0),
    DeltaPrior("r4", 2000.0), LogUniformPrior("ρ4", 1.e-28, 1.e-23), UniformPrior("T4", 0.0, 5.0),
]

prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)
prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]

sampler, results, best_fit_observation = BayesJ.sample(
    observation + background,
    background,
    response_function,
    prior_transform,
    exposure_time,
    exposure_time, # TODO: Make different
    z;
    prior_names=prior_names,
    cluster_model=Model_Piecewise,
    emission_model=emission_model,
    param_wrapper=param_wrapper,
    pixel_edge_angle=pixel_edge_angle,
    centre_radius=centre_radius,
    log_dir="logs/nfw_piecewise",
    # resume="resume",
    integration_limit=integration_limit,
    ultranest_run_args=(
        max_num_improvement_loops=3,
        min_num_live_points=100,)
)

if MPI.Comm_rank(BayesJ.comm) == 0
    best_fit = results["maximum_likelihood"]["point"]
    errlo = results["posterior"]["errlo"]
    errup = results["posterior"]["errup"]

    best_fit_temperature, best_fit_density = Model_Piecewise(param_wrapper(best_fit)[3:end]...)
    mean_fit_temperature, mean_fit_density = Model_Piecewise(param_wrapper(results["posterior"]["mean"])[3:end]...)

    rand_points = [errlo + (errup - errlo) .* rand(Float64, length(errlo)) for i in 1:500]
    models = [Model_Piecewise(param_wrapper(p)[3:end]...) for p in rand_points]
    radii = range(0.0u"kpc", 6000.0u"kpc", length=1000)
    radii_u = ustrip.(u"kpc", radii)

    temperatures = [extrema([model[1](r) for model in models]) for r in radii]
    densities = [extrema([model[2](r) for model in models]) for r in radii]
    output_dir = sampler.logs["run_dir"]

    f = Figure(size=(1200, 1200))
    ax = Axis(f[1, 1], title="Temperature")
    lines!(radii_u, ustrip.(u"keV", temperature.(radii)), label="Original (NFW)")
    lines!(radii_u, ustrip.(u"keV", best_fit_temperature.(radii)), label="Best fit (Piecewise)")
    lines!(radii_u, ustrip.(u"keV", mean_fit_temperature.(radii)), label="Mean fit (Piecewise)")
    band!(radii_u, ustrip.(u"keV", [t[1] for t in temperatures]), ustrip.(u"keV", [t[2] for t in temperatures]), label="68%", alpha=0.5)
    axislegend()
    # ylims!(0.0, 5.0)

    ax2 = Axis(f[2, 1], title="Density")
    lines!(radii_u, ustrip.(u"Msun/kpc^3", density.(radii)), label="Original (NFW)")
    lines!(radii_u, ustrip.(u"Msun/kpc^3", best_fit_density.(radii)), label="Best fit (Piecewise)")
    lines!(radii_u, ustrip.(u"Msun/kpc^3", mean_fit_density.(radii)), label="Mean fit (Piecewise)")
    band!(radii_u, ustrip.(u"Msun/kpc^3", [t[1] for t in densities]), ustrip.(u"Msun/kpc^3", [t[2] for t in densities]), label="68%", alpha=0.5)
    axislegend()

    ax3 = Axis(f[3, 1], title="∝ Pressure")
    lines!(radii_u, ustrip.(u"keV*Msun/kpc^3", temperature.(radii) .* density.(radii)), label="Original (NFW)")
    lines!(radii_u, ustrip.(u"keV*Msun/kpc^3", best_fit_temperature.(radii) .* best_fit_density.(radii)), label="Best fit (Piecewise)")
    lines!(radii_u, ustrip.(u"keV*Msun/kpc^3", mean_fit_temperature.(radii) .* mean_fit_density.(radii)), label="Mean fit (Piecewise)")
    band!(
        radii_u,
        ustrip.(u"keV*Msun/kpc^3", [t[1] * p[1] for (t, p) in zip(temperatures, densities)]),
        ustrip.(u"keV*Msun/kpc^3", [t[2] * p[2] for (t, p) in zip(temperatures, densities)]),
        label="68%", alpha=0.25
    )
    axislegend()

    display(f)


    save("$output_dir/plots/profile.svg", f)
end