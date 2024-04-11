using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using Statistics
include("../src/mpi.jl")
include("../src/params.jl")

# Load Chandra cycle 0 response file
data = FITSData(
    "",
    "",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.arf",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)
response_function, energy_bins, _ = BayesJ.load_response(data, 0.7u"keV", 7.0u"keV")

# Set up emission model
z = 0.1
emission_model = BayesJ.prepare_model_mekal(
    2.e20u"cm^-2",
    energy_bins,
    z,
)

# Set up the full info model
temp_full, density_full = Model_NFW(
    8e15,
    0.13,
    3,
    1.0510,
    5.4905,
    0.3081,
    1.177,
    z=z
)

r_0 = 300u"kpc"
t_0 = uconvert(u"keV", temp_full(r_0))
ρ_0 = uconvert(u"g/cm^3", density_full(r_0))
temp_full, density_full = Model_Constant(r_0 + 100u"kpc", t_0, ρ_0)
println("Temperature: $t_0, Density: $ρ_0")

n_points = 9 # Make sure it's odd so we get r_0 as a point
r_piecewise = collect(LinRange(ustrip(u"kpc", r_0) - 100, ustrip(u"kpc", r_0) + 100, n_points))

# Set up uniform piecewise Model
constant_priors = [DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0)]
for (index, r) in enumerate(r_piecewise)
    global constant_priors
    index -= 1
    if index+1 != Int((n_points+1)/2)
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), DeltaPrior("ρ$index", ustrip(u"g/cm^3", ρ_0)), DeltaPrior("T$index", ustrip(u"keV", t_0))]
    else
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), UniformPrior("ρ$index", 1e-30, 1e-24), UniformPrior("T$index", 1e-4, 1e4)]
    end
    constant_priors = vcat(constant_priors, new_elems)
end

# Set up NFW piecewise model
NFW_priors = [DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0)]
for (index, r) in enumerate(r_piecewise)
    global NFW_priors
    index -= 1
    if index+1 != Int((n_points+1)/2)
        density = density_full(r_piecewise[index+1]u"kpc")
        temp = temp_full(r_piecewise[index+1]u"kpc")
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), 
            DeltaPrior("ρ$index", ustrip(u"g/cm^3", density)), 
            DeltaPrior("T$index", ustrip(u"keV", temp))]
    else
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), UniformPrior("ρ$index", 1e-30, 1e-24), UniformPrior("T$index", 1e-4, 1e4)]
    end
    NFW_priors = vcat(NFW_priors, new_elems)
end


# For convenience we want the pixel width to be a factor of 100kpc. These params give ~24.95kpc, so center radius as 4 will block the 100kpc range.
shape = (32, 32)
pixel_edge_angle = 13.1u"arcsecondᵃ"
exposure_time = 5.0e6u"s"
centre_radius = 4
integration_limit = 10u"Mpc"

pixel_size = ustrip(u"kpc", angle_to_length(pixel_edge_angle, z)) # Approx 24.95

# Now we can make an observation with the full model to see what it looks like
function observations_from_profiles(t, d)
    predicted_count_rate = BayesJ.make_observation(
            temp_full,
            density_full,
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
    replace!(predicted_count_rate, missing=>0)

    # Random.seed!(42)
    observation = pois_rand.(predicted_count_rate)
    bg_rate = 0.01*maximum(predicted_count_rate)
    println("Background Rate $bg_rate")
    background = Array{Int64}(undef, size(observation)...)
    for i in eachindex(background)
        background[i] = pois_rand(bg_rate) + 1
    end
    observation += background
    noise = observation - predicted_count_rate
    return (predicted_count_rate, observation, bg_rate, noise)
end

expected_rate, observed_rate, background_rate, noise_rate = observations_from_profiles(temp_full, density_full)

# Plot temp, density profiles
function plot_profiles(do_save=false)
    # Plotting variables
    radii = 1:500
    t(r) = ustrip(u"keV", temp_full(r * 1u"kpc"))
    d(r) = ustrip(u"g/cm^3", density_full(r * 1u"kpc"))
    t_piecewise(r_piecewise) = ustrip(u"keV", temp_full(r_piecewise * 1u"kpc"))
    d_piecewise(r_piecewise) = ustrip(u"g/cm^3", density_full(r_piecewise * 1u"kpc"))
    
    f = Figure()
    # Temperature plot
    ax = Axis(f[1, 1], ylabel="Temperature [keV]")
    hidexdecorations!(ax, grid=false)
    lines!(radii, t.(radii), label="Full Model")
    scatter!(r_piecewise, fill(ustrip(u"keV", t_0), n_points), label="Uniform Model")
    scatter!(r_piecewise, t_piecewise, label="NFW Model")
    vlines!(ax, ustrip(u"kpc", r_0), color=:black, label=L"r_0")

    # Density plot
    Label(f[2, 1, TopLeft()], "Density", font=:bold, padding = (0, 0, -30, 0))
    ax2 = Axis(f[2, 1], xlabel="Radius [kpc]", ylabel="Density [g/cm^3]")
    lines!(radii, d.(radii))
    scatter!(r_piecewise, fill(ustrip(u"g/cm^3", ρ_0), n_points))
    scatter!(r_piecewise, d_piecewise)
    vlines!(ax2, ustrip(u"kpc", r_0), color=:black)

    # Plot details 
    Label(f[1, 1, Top()], "Test Models Temperature and Density profiles", font=:bold, fontsize=20, padding = (0, 0, 20, 0))
    f[1:2, 2] = Legend(f, ax, framevisible=false)
    
    if do_save
        display(f)
        save("../signal_noise_tests/model_profiles.png", f)
    end
end

# Plot cluster counts
function plot_cluster_counts(do_save=false)
    f = Figure(size=(800,800))
    x = [pixel_size * (pix - (shape[1]+1)/2) for pix in 1:shape[1]]
    y = [pixel_size * (pix - (shape[2]+1)/2) for pix in 1:shape[2]]
    expected_counts = [sum(expected_rate[:, Int(x), Int(y)]) for x in 1:shape[1], y in 1:shape[2]]
    observed_counts = [sum(observed_rate[:, Int(x), Int(y)]) for x in 1:shape[1], y in 1:shape[2]]
    noise_counts = expected_counts - observed_counts

    ax = Axis(f[1, 1], title="Expected Counts", ylabel="Distance from center [kpc]", yticks=-400:50:400)
    hm = heatmap!(ax, x, y, expected_counts)
    Colorbar(f[1, 2], hm)

    ax2 = Axis(f[1, 3], title="Sample Observation Counts")
    hm = heatmap!(ax2, x, y, observed_counts)
    Colorbar(f[1, 4], hm)

    ax3 = Axis(f[2, 1], title="Noise Counts", xlabel="Distance from center [kpc]", xticks=-400:50:400, xticklabelrotation=pi/2,
        ylabel="Distance from center [kpc]", yticks=-400:50:400)
    hm = heatmap!(ax3, x, y, noise_counts)
    Colorbar(f[2, 2], hm)

    ax4 = Axis(f[2, 3], title="Signal/Noise", xlabel="Distance from center [kpc]", xticks=-400:50:400, xticklabelrotation=pi/2)
    hm = heatmap!(ax4, x, y, observed_counts./noise_counts)
    Colorbar(f[2, 4], hm)
    hidexdecorations!(ax)
    hidedecorations!(ax2)
    hideydecorations!(ax4)
    if do_save
        display(f)
        save("../signal_noise_tests/cluster_counts.svg", f)
    end
end 

# Plot spectrumn averaged within some annulus
function plot_annulus_spectra(do_save=false)
    x = [pixel_size * (pix - (shape[1]+1)/2) for pix in 1:shape[1]]
    y = [pixel_size * (pix - (shape[2]+1)/2) for pix in 1:shape[2]]
    annulus_pixel_mask = similar(expected_rate, Bool)
    for i in 1:shape[1], j in 1:shape[2]
        x_val = x[i]
        y_val = x[j]
        annulus_pixel_mask[:, i, j] .= hypot(x_val, y_val)>200 && hypot(x_val, y_val)<400 ? true : false
    end
    function annulus_spectrum(observation, mask)
        annulus_pixels = reshape(observation[mask], size(observation)[1], :)
        mean(annulus_pixels, dims=2)[:]
    end
    temp_uniform, density_uniform = Model_Constant(400u"kpc", t_0, ρ_0)
    expected_uniform, observed_uniform, background_uniform, noise_uniform = observations_from_profiles(temp_uniform, density_uniform)
    f = Figure(size=(800,800))
    # ax = Axis(f[1, 1], title="Annulus", ylabel="Distance from center [kpc]", xlabel="Distance from center [kpc]", xticks=-400:50:400, 
    #     yticks=-400:50:400, xticklabelrotation=pi/2)
    # hm = heatmap!(ax, x, y, annulus_pixel_mask[1, :, :])
    # Colorbar(f[1, 2], hm)

    ax = Axis(f[1, 1], title="Averaged observation annulus spectrum", xlabel="Channel", ylabel="Counts")
    lines!(ax, annulus_spectrum(expected_uniform, annulus_pixel_mask), label="Expected Uniform")
    lines!(ax, annulus_spectrum(observed_uniform, annulus_pixel_mask), label="Observation Uniform")
    # lines!(ax, annulus_spectrum(expected_rate, annulus_pixel_mask), label="Expected NFW")
    # lines!(ax, annulus_spectrum(observed_rate, annulus_pixel_mask), label="Observation NFW")
    axislegend(ax)

    ax = Axis(f[1, 2], title="Averaged signal/noise annulus spectrum", xlabel="Channel", ylabel="Ratio")
    lines!(ax, annulus_spectrum(observed_uniform, annulus_pixel_mask)./annulus_spectrum(noise_uniform, annulus_pixel_mask), label="Uniform")
    # lines!(ax, annulus_spectrum(observed_rate, annulus_pixel_mask)./annulus_spectrum(noise_rate, annulus_pixel_mask), label="NFW")
    axislegend(ax)

    if do_save
        display(f)
        save("../signal_noise_tests/annulus_spectra.png", f)
    end
end 

plot_profiles(true)
plot_cluster_counts(true)
plot_annulus_spectra(true)

# temp plot function 
function plot_res(obs::Array, pred::Array, path::AbstractString)
    if BayesJ.isroot()
        if pred[:, 1, 1] != pred[:, 4, 4]
            @warn "Something funky with the best fit array"
        end
        slice_obs = Vector{Float64}(obs[:, 4, 4])
        slice_pred = Vector{Float64}(pred[:, 4, 4])

        f = Figure()
        ax = Axis(f[1, 1], xlabel="Channel", ylabel="Counts", yscale=Makie.pseudolog10)
        lines!(slice_obs ./ maximum(slice_obs), label="Observation")
        lines!(slice_pred ./ maximum(slice_pred), label="Best Fit")
        axislegend()

        save(joinpath(path, "spectra.svg"), f)
    end
end

# Run model fitting and export
function run_model_fit(obs, bg, priors, log_dir=nothing)
    prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)
    prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]

    sampler, result, best_fit_observation = BayesJ.sample(
        obs,
        bg,
        response_function,
        prior_transform,
        exposure_time,
        exposure_time, 
        z;
        prior_names=prior_names,
        cluster_model=Model_Piecewise,
        emission_model=emission_model,
        param_wrapper=param_wrapper,
        pixel_edge_angle=pixel_edge_angle,
        centre_radius=centre_radius,
        log_dir=log_dir,
        integration_limit=integration_limit,
        ultranest_run_args=(
            max_num_improvement_loops=3,
            min_num_live_points=100,
            show_status=false,
            viz_callback=false
        )
    )

    try
        output_dir = sampler.logs["run_dir"]
        if output_dir isa AbstractString
            plot_res(obs, best_fit_observation, joinpath(output_dir, "plots"))
        end
    catch e
        if !(e isa KeyError)
            rethrow(e)
        end
        BayesJ.@mpidebug "Skipping spectra plot because unable to find log information"
    end
end

# for i in 1:10
#     global expected_rate, observed_rate, background_rate, noise_rate
#     if isroot()
#         expected_rate, observed_rate, background_rate, noise_rate = observations_from_profiles(temp_full, density_full)
#     end
#     observed_rate = MPI.bcast(observed_rate, 0, comm)
#     background_rate = MPI.bcast(background_rate, 0, comm)
#     print("Rank: $(MPI.Comm_rank(comm)) has $(background_rate[1:10, 1, 1])")
#     run_model_fit(observed_rate, background_rate, NFW_priors, "../signal_noise_tests/NFW_annulus")
#     observed_rate = MPI.bcast(observed_rate, 0, comm)
#     background_rate = MPI.bcast(background_rate, 0, comm)
#     run_model_fit(observed_rate, background_rate, constant_priors, "../signal_noise_tests/constant_annulus")
# end