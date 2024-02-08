using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using MPI
using LinearAlgebra: I

z = 0.5
shape = (32, 32)

data = FITSData(
    "",
    "",
    "data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.arf",
    "data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

response_function, energy_bins, _ = BayesJ.load_response(data, 0.7u"keV", 7.0u"keV")

pixel_edge_angle = 4.92u"arcsecondᵃ"
exposure_time = 30.0e6u"s"
centre_radius = 0

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
)

@assert all(isfinite, predicted_count_rate)

observation = pois_rand.(predicted_count_rate)
background_rate = rand() * min(1, maximum(predicted_count_rate))
background = Array{Int64}(undef, size(observation)...)
for i in eachindex(background)
    background[i] = pois_rand(background_rate) + 1
end

@assert all(isfinite, observation)
@assert all(isfinite, background)

if BayesJ.isroot()
    display(lines(predicted_count_rate[:, 24, 24]))
    display(lines(observation[:, 24, 24]))
end

priors = [
    DeltaPrior("x0", 0.0),
    DeltaPrior("y0", 0.0),
    UniformPrior("MT_500", 1.0e14, 4.0e14),
    UniformPrior("fg_500", 0.08, 0.5),
    DeltaPrior("c_500", 3.0),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
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
    cluster_model=Model_NFW,
    emission_model=emission_model,
    param_wrapper=param_wrapper,
    pixel_edge_angle=pixel_edge_angle,
    centre_radius=centre_radius,
    log_dir="logs/nfw_piecewise",
    # resume="resume",
    ultranest_run_args=(
        max_num_improvement_loops=3,
        min_num_live_points=100,)
)

if MPI.Comm_rank(BayesJ.comm) == 0
    best_fit = results["maximum_likelihood"]["point"]
    errlo = results["posterior"]["errlo"]
    errup = results["posterior"]["errup"]

    best_fit_temperature, best_fit_density = Model_NFW(param_wrapper(best_fit)[3:end]..., z=z)
    mean_fit_temperature, mean_fit_density = Model_NFW(param_wrapper(results["posterior"]["mean"])[3:end]..., z=z)

    rand_points = [errlo + (errup - errlo) .* rand(Float64, length(errlo)) for i in 1:500]
    models = [Model_Piecewise(param_wrapper(p)[3:end]...) for p in rand_points]
    radii = range(0.0u"kpc", 1000u"kpc", length=1000)
    radii_u = ustrip.(u"kpc", radii)

    temperatures = [extrema([model[1](r) for model in models]) for r in radii]
    densities = [extrema([model[2](r) for model in models]) for r in radii]
    output_dir = sampler.logs["run_dir"]

    f = Figure(size=(1200, 1200))
    ax = Axis(f[1, 1], title="Temperature", xlabel="Radius (kpc)", ylabel="Temperature (keV)")
    lines!(radii_u, ustrip.(u"keV", temperature.(radii)), label="Original (NFW)")
    lines!(radii_u, ustrip.(u"keV", best_fit_temperature.(radii)), label="Best fit (Piecewise)")
    lines!(radii_u, ustrip.(u"keV", mean_fit_temperature.(radii)), label="Mean fit (Piecewise)")
    band!(radii_u, ustrip.(u"keV", [t[1] for t in temperatures]), ustrip.(u"keV", [t[2] for t in temperatures]), label="68%", alpha=0.5)
    axislegend()
    # ylims!(0.0, 5.0)

    ax2 = Axis(f[2, 1], title="Density", xlabel="Radius (kpc)", ylabel="Density (Msun/kpc^3)", yscale=Makie.pseudolog10)
    lines!(radii_u, ustrip.(u"Msun/kpc^3", density.(radii)), label="Original (NFW)")
    lines!(radii_u, ustrip.(u"Msun/kpc^3", best_fit_density.(radii)), label="Best fit (Piecewise)")
    lines!(radii_u, ustrip.(u"Msun/kpc^3", mean_fit_density.(radii)), label="Mean fit (Piecewise)")
    band!(radii_u, ustrip.(u"Msun/kpc^3", [t[1] for t in densities]), ustrip.(u"Msun/kpc^3", [t[2] for t in densities]), label="68%", alpha=0.5)
    axislegend()

    ax3 = Axis(f[3, 1], title="∝ Pressure", xlabel="Radius (kpc)", ylabel="∝ Pressure (kev Msun / kpc^3)", yscale=Makie.pseudolog10)
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