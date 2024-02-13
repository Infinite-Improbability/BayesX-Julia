using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using LinearAlgebra: I

z = 0.1
shape = (64, 64)

# data = FITSData(
#     "",
#     "",
#     "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.arf",
#     "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.rmf",
#     0.492u"arcsecondᵃ"
# )

# response_function, energy_bins, _ = BayesJ.load_response(data, 0.7u"keV", 7.0u"keV")

energy_bins = range(0.7u"keV", 7.0u"keV", length=3)
response_function = 250u"cm^2" * Matrix(I, length(energy_bins) - 1, length(energy_bins) - 1)

pixel_edge_angle = 20.0u"arcsecondᵃ"
exposure_time = 3.0e9u"s"
centre_radius = 0
integration_limit = 875u"kpc"

emission_model = BayesJ.prepare_model_mekal(
    0.0u"cm^-2",
    energy_bins,
    z,
)

cluster_model(args...; kwargs...) = Model_Piecewise(args...; kwargs...)

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
    limit=integration_limit
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
@assert !all(iszero, observation) "Maximum pcr is $(maximum(predicted_count_rate))"

t(r) = ustrip(u"keV", temperature(r * 1u"kpc"))
d(r) = ustrip(u"g/cm^3", density(r * 1u"kpc"))

r1 = 1000
r2 = 1100
r3 = 1200
r4 = 1250
r5 = 1300
r6 = 1350
r7 = 1400
r8 = 1450
r9 = 1500

@mpiinfo "Target values" r2 t(r2) d(r2)

if BayesJ.isroot()
    pixel_edge_length = BayesJ.angle_to_length(pixel_edge_angle, z)

    halfx = size(observation, 2) ÷ 2
    halfy = size(observation, 3) ÷ 2
    offset_x = r2 / ustrip(u"kpc", pixel_edge_length) + halfx
    offset_y = r2 / ustrip(u"kpc", pixel_edge_length) + halfy
    rx = round(Int, offset_x)
    ry = round(Int, offset_y)

    flux = Vector{Cfloat}(undef, size(response_function, 2))
    emission_model(flux, temperature(r2 * 1u"kpc"), density(r2 * 1u"kpc"))
    @mpiinfo "Flux estimate" extrema(flux) flux
    @assert sum(observation[:, rx, ry]) > 0 "maximum pcr is $(maximum(predicted_count_rate[:, rx, ry]))"

    f = Figure()
    ax = Axis(f[1, 1], yscale=Makie.pseudolog10)
    lines!(observation[:, rx, ry], color=:green)
    display(f)
end

priors = [
    DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0),
    DeltaPrior("r10", 10.0), DeltaPrior("ρ10", d(10.0)), DeltaPrior("T10", t(10.0)),
    DeltaPrior("r50", 50.0), DeltaPrior("ρ50", d(50.0)), DeltaPrior("T50", t(50.0)),
    DeltaPrior("r100", 100.0), DeltaPrior("ρ100", d(100.0)), DeltaPrior("T100", t(100.0)),
    DeltaPrior("r200", 200.0), DeltaPrior("ρ200", d(200.0)), DeltaPrior("T200", t(200.0)),
    DeltaPrior("r300", 300.0), DeltaPrior("ρ300", d(300.0)), DeltaPrior("T300", t(300.0)),
    DeltaPrior("r400", 400.0), DeltaPrior("ρ400", d(400.0)), DeltaPrior("T400", t(400.0)),
    DeltaPrior("r500", 500.0), DeltaPrior("ρ500", d(500.0)), DeltaPrior("T500", t(500.0)),
    DeltaPrior("r550", 550.0), DeltaPrior("ρ550", d(550.0)), DeltaPrior("T550", t(550.0)),
    DeltaPrior("r600", 600.0), DeltaPrior("ρ600", d(600.0)), DeltaPrior("T600", t(600.0)),
    DeltaPrior("r650", 650.0), DeltaPrior("ρ650", d(650.0)), DeltaPrior("T650", t(650.0)),
    DeltaPrior("r700", 700.0), DeltaPrior("ρ700", d(700.0)), DeltaPrior("T700", t(700.0)),
    DeltaPrior("r750", 750.0), DeltaPrior("ρ750", d(750.0)), DeltaPrior("T750", t(750.0)),
    DeltaPrior("r800", 800.0), DeltaPrior("ρ800", d(800.0)), DeltaPrior("T800", t(800.0)),
    DeltaPrior("r850", 850.0), DeltaPrior("ρ850", d(850.0)), DeltaPrior("T850", t(850.0)),
    DeltaPrior("r900", 900.0), DeltaPrior("ρ900", d(900.0)), DeltaPrior("T900", t(900.0)),
    DeltaPrior("r950", 950.0), DeltaPrior("ρ950", d(950.0)), DeltaPrior("T950", t(950.0)),
    DeltaPrior("r1", r1), DeltaPrior("ρ1", d(r1)), DeltaPrior("T1", t(r1)),
    DeltaPrior("r2", r2), DeltaPrior("ρ2", d(r2)), LogUniformPrior("T2", 0.5 * t(r3), 2 * t(r1)),
    DeltaPrior("r3", r3), DeltaPrior("ρ3", d(r3)), DeltaPrior("T3", t(r3)),
    DeltaPrior("r4", r4), DeltaPrior("ρ4", d(r4)), DeltaPrior("T4", t(r4)),
    DeltaPrior("r5", r5), DeltaPrior("ρ5", d(r5)), DeltaPrior("T5", t(r5)),
    DeltaPrior("r6", r6), DeltaPrior("ρ6", d(r6)), DeltaPrior("T6", t(r6)),
    DeltaPrior("r7", r7), DeltaPrior("ρ7", d(r7)), DeltaPrior("T7", t(r7)),
    DeltaPrior("r8", r8), DeltaPrior("ρ8", d(r8)), DeltaPrior("T8", t(r8)),
    DeltaPrior("r9", r9), DeltaPrior("ρ9", d(r9)), DeltaPrior("T9", t(r9)),
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
    cluster_model=cluster_model,
    emission_model=emission_model,
    param_wrapper=param_wrapper,
    pixel_edge_angle=pixel_edge_angle,
    centre_radius=centre_radius,
    log_dir="../logs/pw_piecewise/",
    integration_limit=integration_limit,
    # resume="resume",
    ultranest_run_args=(
        max_num_improvement_loops=1,
        min_num_live_points=100,)
)

@mpiinfo "Making plots"

if BayesJ.isroot()
    best_fit = results["maximum_likelihood"]["point"]
    errlo = results["posterior"]["errlo"]
    errup = results["posterior"]["errup"]

    best_fit_temperature, best_fit_density = cluster_model(param_wrapper(best_fit)[3:end]..., z=z)
    mean_fit_temperature, mean_fit_density = cluster_model(param_wrapper(results["posterior"]["mean"])[3:end]..., z=z)

    rand_points = [errlo + (errup - errlo) .* rand(Float64, length(errlo)) for i in 1:500]
    models = [cluster_model(param_wrapper(p)[3:end]..., z=z) for p in rand_points]
    radii = range(0.0u"kpc", r9 * 1u"kpc", length=1000)
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
    axislegend(position=:lb)
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