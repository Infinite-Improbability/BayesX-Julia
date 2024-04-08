using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using Statistics
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
    8e13,
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
println("Temperature: $t_0, Density: $ρ_0")

n_points = 9 # Make sure it's odd so we get r_0 as a point
r_piecewise = collect(LinRange(ustrip(u"kpc", r_0) - 100, ustrip(u"kpc", r_0) + 100, n_points))

# Set up uniform piecewise Model
constant_priors = [DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0)]
for (index, r) in enumerate(r_piecewise)
    global constant_priors
    index -= 1
    if index+1 != Int((n_points+1)/2)
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), DeltaPrior("ρ$index", ρ_0), DeltaPrior("T$index", t_0)]
    else
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), UniformPrior("ρ$index", 1e-40, 1e-34), UniformPrior("T$index", 1e-10, 1e-4),]
    end
    constant_priors = vcat(constant_priors, new_elems)
end

# Set up NFW piecewise model
NFW_priors = [DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0)]
for (index, r) in enumerate(r_piecewise)
    global NFW_priors
    index -= 1
    if index+1 != Int((n_points+1)/2)
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), 
            DeltaPrior("ρ$index", density_full(r_piecewise[index+1]u"kpc")), 
            DeltaPrior("T$index", temp_full(r_piecewise[index+1]u"kpc"))]
    else
        new_elems = [DeltaPrior("r$index", r_piecewise[index+1]), UniformPrior("ρ$index", 1e-40, 1e-34), UniformPrior("T$index", 1e-10, 1e-4),]
    end
    NFW_priors = vcat(NFW_priors, new_elems)
end

# Plotting
radii = 1:500
t(r) = ustrip(u"keV", temp_full(r * 1u"kpc"))
d(r) = ustrip(u"g/cm^3", density_full(r * 1u"kpc"))
t_piecewise(r_piecewise) = ustrip(u"keV", temp_full(r_piecewise * 1u"kpc"))
d_piecewise(r_piecewise) = ustrip(u"g/cm^3", density_full(r_piecewise * 1u"kpc"))

f = Figure()
# Temperature plot
# Label(f[1, 1, TopLeft()], "Temperature", font=:bold, padding = (0, 0, -30, 0))
ax = Axis(f[1, 1], ylabel="Temperature [keV]")
hidexdecorations!(ax, grid=false)
lines!(radii, t.(radii), label="Full Model")
scatter!(r_piecewise, fill(ustrip(u"keV", t_0), n_points), label="Uniform Model")
scatter!(r_piecewise, t_piecewise, label="NFW Model")
vlines!(ax, ustrip(u"kpc", r_0), color=:black, label=L"r_0")

# Density plot
Label(f[2, 1, TopLeft()], "Density", font=:bold, padding = (0, 0, -30, 0))
ax2 = Axis(f[2, 1], xlabel="Radius [kpc]", ylabel="Density [g/cm^3]", yscale=log10)
lines!(radii, d.(radii))
scatter!(r_piecewise, fill(ustrip(u"g/cm^3", ρ_0), n_points))
scatter!(r_piecewise, d_piecewise)
vlines!(ax2, ustrip(u"kpc", r_0), color=:black)

# Plot details 
Label(f[1, 1, Top()], "Test Models Temperature and Density profiles", font=:bold, fontsize=20, padding = (0, 0, 20, 0))
f[1:2, 2] = Legend(f, ax, framevisible=false)
# display(f)
# save("../signal_noise_tests/model_profiles.svg", f)

# For convenience we want the pixel width to be a factor of 100kpc. These params give ~24.95kpc, so center radius as 4 will block the 100kpc range.
shape = (32, 32)
pixel_edge_angle = 13.1u"arcsecondᵃ"
exposure_time = 3.0e8u"s"
centre_radius = 4
integration_limit = 10u"Mpc"

pixel_size = ustrip(u"kpc", angle_to_length(pixel_edge_angle, z)) # Approx 24.95
# println("Pixel size [kpc]: $pixel_size")
# Now we can make an observation with the full model to see what it looks like
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

Random.seed!(42)
observation = pois_rand.(predicted_count_rate)
bg_rate = rand() * min(1, maximum(predicted_count_rate))
background = Array{Int64}(undef, size(observation)...)
for i in eachindex(background)
    background[i] = pois_rand(bg_rate) + 1
end
observation += background
noise = observation - predicted_count_rate

# Plot cluster counts
f = Figure(size=(800,800))
x = [pixel_size * (pix - (shape[1]+1)/2) for pix in 1:shape[1]]
y = [pixel_size * (pix - (shape[2]+1)/2) for pix in 1:shape[2]]
expected_counts = [sum(predicted_count_rate[:, Int(x), Int(y)]) for x in 1:shape[1], y in 1:shape[2]]
observed_counts = [sum(observation[:, Int(x), Int(y)]) for x in 1:shape[1], y in 1:shape[2]]
noise_counts = expected_counts - observed_counts

ax = Axis(f[1, 1], title="Expected Counts", ylabel="Distance from center [kpc]", yticks=[round(y_val) for y_val in y[1:2:end]])
hm = heatmap!(ax, x, y, expected_counts)
Colorbar(f[1, 2], hm)

ax2 = Axis(f[1, 3], title="Sample Observation Counts")
hm = heatmap!(ax2, x, y, observed_counts)
Colorbar(f[1, 4], hm)

ax3 = Axis(f[2, 1], title="Noise Counts", xlabel="Distance from center [kpc]", xticks=[round(x_val) for x_val in x[1:2:end]], xticklabelrotation=pi/2,
    ylabel="Distance from center [kpc]", yticks=[round(y_val) for y_val in y[1:2:end]])
hm = heatmap!(ax3, x, y, noise_counts)
Colorbar(f[2, 2], hm)

ax4 = Axis(f[2, 3], title="Signal/Noise", xlabel="Distance from center [kpc]", xticks=[round(x_val) for x_val in x[1:2:end]], xticklabelrotation=pi/2)
hm = heatmap!(ax4, x, y, observed_counts./noise_counts)
Colorbar(f[2, 4], hm)
hidexdecorations!(ax)
hidedecorations!(ax2)
hideydecorations!(ax4)
# display(f)
# save("../signal_noise_tests/cluster_counts.svg", f)

# Plot spectrumn averaged within some annulus
annulus_pixel_mask = similar(predicted_count_rate, Bool)
for i in 1:shape[1], j in 1:shape[2]
    x_val = x[i]
    y_val = x[j]
    annulus_pixel_mask[:, i, j] .= hypot(x_val, y_val)>200 && hypot(x_val, y_val)<400 ? true : false
end

function annulus_spectrum(observation, mask)
    annulus_pixels = reshape(observation[mask], size(observation)[1], :)
    mean(annulus_pixels, dims=2)[:]
end

f = Figure(size=(800,800))

ax = Axis(f[1, 1], title="Selected area for averaging", ylabel="Distance from center [kpc]", yticks=[round(y_val) for y_val in y[1:2:end]])
hm = heatmap!(ax, x, y, annulus_pixel_mask[1, :, :])
Colorbar(f[1, 2], hm)

ax = Axis(f[2, 1], title="Averaged annulus spectrum")
lines!(ax, annulus_spectrum(observation, annulus_pixel_mask))
display(f)
save("../signal_noise_tests/annulus_spectra.svg", f)
