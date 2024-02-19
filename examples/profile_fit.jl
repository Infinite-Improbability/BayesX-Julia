using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using LinearAlgebra: I

z = 0.1
shape = (32, 32)

data = FITSData(
    "",
    "",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.arf",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

response_function, energy_bins, _ = BayesJ.load_response(data, 0.7u"keV", 7.0u"keV")

# energy_bins = range(0.7u"keV", 7.0u"keV", step=0.01u"keV")
# response_function = 250u"cm^2" * Matrix(I, length(energy_bins) - 1, length(energy_bins) - 1)

pixel_edge_angle = 20.0u"arcsecondᵃ"
exposure_time = 3.0e8u"s"
centre_radius = 1
integration_limit = 10u"Mpc"

emission_model = BayesJ.prepare_model_mekal(
    0.0u"cm^-2",
    energy_bins,
    z,
)

cluster_model(args...; kwargs...) = Model_Piecewise(args...; kwargs...)

temperature, density = Model_Vikhlinin2006(
    4.705e-3,
    0.247e-1,
    94.6,
    75.83,
    0.916,
    0.526,
    3.607,
    4.943,
    1239.9,
    3.61,
    0.27,
    57.0,
    3.88,
    1420.0,
    0.12,
    5.0,
    10.0,
    0.0,
    z=z
)

# radii = 1:2000
# t(r) = ustrip(u"keV", temperature(r * 1u"kpc"))
# d(r) = ustrip(u"g/cm^3", density(r * 1u"kpc"))

# f = Figure()
# ax = Axis(f[1, 1], title="Temperature", xlabel="Radius (kpc)", ylabel="Temperature (keV)")
# lines!(radii, t.(radii))
# ax2 = Axis(f[2, 1], title="Density", xlabel="Radius (kpc)", ylabel="Density (g/cm^3)", yscale=log10)
# lines!(radii, d.(radii))
# display(f)

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

replace!(predicted_count_rate, missing => 0.0)
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

t(r) = ustrip(u"keV", temperature(r))
d(r) = ustrip(u"g/cm^3", density(r))

r1 = 500u"kpc"
r2 = 600u"kpc"
r3 = 700u"kpc"

@mpiinfo "Target values" r2 t(r2) d(r2)

if BayesJ.isroot()
    pixel_edge_length = BayesJ.angle_to_length(pixel_edge_angle, z)

    halfx = size(observation, 2) ÷ 2
    halfy = size(observation, 3) ÷ 2
    offset_x = ustrip(u"kpc/kpc", r2 / pixel_edge_length) + halfx
    offset_y = ustrip(u"kpc/kpc", r2 / pixel_edge_length) + halfy
    rx = round(Int, offset_x)
    ry = round(Int, offset_y)

    flux = Vector{Cfloat}(undef, size(response_function, 2))
    emission_model(flux, temperature(r2), density(r2))
    @mpiinfo "Flux estimate" extrema(flux) flux
    @assert sum(observation[:, rx, ry]) > 0 "maximum pcr is $(maximum(predicted_count_rate[:, rx, ry]))"

    f = Figure()
    ax = Axis(f[1, 1], yscale=Makie.pseudolog10)
    lines!(observation[:, rx, ry], color=:green)
    display(f)
end

priors = [
    DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0),
    DeltaPrior("r2", ustrip(u"kpc", r2)), UniformPrior("ρ2", 0.5 * d(r3), 2 * d(r1)), UniformPrior("T2", 0.5 * t(r3), 2 * t(r1)),
]

prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)
prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]

function cluster_model(ri, ρ, T; kwargs...)
    let temperature = temperature, density = density, r1 = r1, r3 = r3
        temp, den = Model_Piecewise(ustrip(u"kpc", r1), d(r1), t(r1), ri, ρ, T, ustrip(u"kpc", r3), d(r3), t(r3))

        function custom_temperature(r)
            let r1 = r1, r3 = r3, temperature = temperature
                if r1 < r < r3
                    return temp(r)
                else
                    return temperature(r)
                end
            end
        end

        function custom_density(r)
            let r1 = r1, r3 = r3, density = density
                if r1 < r < r3
                    return den(r)
                else
                    return density(r)
                end
            end
        end

        return custom_temperature, custom_density
    end
end

sampler, results, _ = BayesJ.sample(
    observation + background,
    background,
    response_function,
    prior_transform,
    exposure_time,
    exposure_time,
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
    radii = range(0.01u"kpc", 1.2 * r2, length=1000)
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

    ax2 = Axis(f[2, 1], title="Density", xlabel="Radius (kpc)", ylabel="Density (Msun/kpc^3)", yscale=Makie.pseudolog10) \
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