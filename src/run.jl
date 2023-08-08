using Unitful, DimensionfulAngles
using Plots
using Profile

using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")
stepsampler = pyimport_conda("ultranest.stepsampler", "ultranest")

export sample

include("mpi.jl")
include("gas_models.jl")
include("io.jl")
include("likelihood.jl")

"""
    sample(observed::Array, observed_background::Array, response_function::Matrix, priors::PriorSet, obs_exposure_time, bg_exposure_time, redshift; emission_model, pixel_edge_angle, background_rate, average_effective_area, center_radius)

Configure some necessary variables and launch ultranest.

* The observed array includes the background.
* The response function includes both the RMF and ARF, as described in `apply_response_function`.
* The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
* Prior names should match arguments for the gas model.
* The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
 This area is assumed to be square with the edge angle giving the side length.
* The average effective area is the effective area of the telescope averaged across energies,
 used with the total background rate across all channels (counts per unit telescope area per sky angle per second)
 to calculate the background counts per second per channel per pixel.
 * `center_radius` controls how many pixels are excluded from the core.
"""
function sample(
    observed::T,
    observed_background::T,
    response_function::Matrix,
    priors::PriorSet,
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
    redshift::Real;
    emission_model,
    pixel_edge_angle=0.492u"arcsecondᵃ",
    background_rate=8.4e-6u"cm^-2/arcminuteᵃ^2/s",
    average_effective_area=250u"cm^2",
    center_radius
) where {T<:AbstractArray}
    @mpidebug "Preparing for ultranest"

    predicted_bg_rate = background_rate / size(observed)[1] * average_effective_area * pixel_edge_angle^2
    predicted_obs_bg = predicted_bg_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = predicted_bg_rate * bg_exposure_time # Used for log likelihood

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)

    @assert all(isfinite, observed)
    @assert all(isfinite, observed_background)
    @assert all(isfinite, log_obs_factorial)

    @assert size(observed) == size(observed_background)

    shape = [i for i in size(observed)][2:3]

    @mpiinfo "Observation has shape $(size(observed))"
    @mpiinfo "Background has shape $(size(observed_background))"
    @mpiinfo "Response matrix has shape $(size(response_function))"

    # a wrapper to handle running the gas model and likelihood calculation
    @mpidebug "Generating likelihood wrapper"
    function likelihood_wrapper(params)
        @mpirankeddebug "Likelihood wrapper called" params

        predicted = Model_NFW_GNFW(
            params[1],
            params[2],
            1.0510, # Using universal values from Arnaud 2010
            5.4905,
            0.3081,
            1.177,
            redshift,
            shape,
            pixel_edge_angle,
            emission_model,
            obs_exposure_time,
            response_function,
            # (params[3], params[4])
            (0, 0),
            center_radius
        ) .+ predicted_obs_bg


        @mpirankeddebug "Predicted results generated"

        # @mpirankedinfo "Taking heap snapshot"
        # Profile.take_heap_snapshot()
        # @mpirankedinfo "Snapshot made"
        # [display(heatmap(dropdims(sum(p, dims=1), dims=1))) for p in predicted]

        return log_likelihood(
            observed,
            observed_background,
            predicted,
            predicted_bg_bg,
            log_obs_factorial
        )
    end

    # ultranest setup
    @mpidebug "Creating sampler"
    # paramnames = ["MT_200", "fg_200", "x_0", "y_0"] # move to pairs with prior objects?
    paramnames = ["MT_200", "fg_200"]
    sampler = ultranest.ReactiveNestedSampler(
        paramnames,
        likelihood_wrapper,
        transform=priors.transform,
        vectorized=false,
        log_dir="logs"
    )

    # @mpidebug "Creating stepsampler"
    # sampler.stepsampler = stepsampler.SliceSampler(
    #     nsteps=25,
    #     generate_direction=stepsampler.generate_mixture_random_direction,
    # )

    # run Ultranest
    @mpiinfo "Launching sampler"
    results = sampler.run()

    # output data
    @mpidebug "Sampler done"
    # print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot_corner()
    sampler.plot_trace()
    sampler.plot_run()

    return (sampler, results)
end

"""
    sample(data::Dataset, energy_range, priors::Dict, nhCol, redshift; bin_size, use_interpolation, center_radius)

Run Bayesian inference on a given set of `data`, considering only the selected energy range.

* `bin_size` controls spatial size of bins, in pixels
* `use_interpolation` controls whether the gas emission model uses interpolation
* `center_radius` controls the excluded area in the core, in pixels
"""
function sample(
    data::Dataset,
    energy_range::AbstractRange{<:Unitful.Energy},
    priors::Dict{<:AbstractString,<:Prior},
    nHcol::SurfaceDensity,
    redshift::Real;
    bin_size::Real=10,
    use_interpolation::Bool=true,
    center_radius=4
)
    @mpiinfo "Loading data"

    observation, observed_background = load_data(data)

    obs = bin_events(data, observation.first, energy_range, 3700:bin_size:4200, 4100:bin_size:4550)
    bg = bin_events(data, observed_background.first, energy_range, 3700:bin_size:4200, 4100:bin_size:4550)
    pixel_edge_angle = bin_size * data.pixel_edge_angle
    @mpidebug "Done binning events"

    @assert size(obs) == size(bg)

    @mpidebug "Making transform"
    prior_set = generate_transform(priors)

    @mpidebug "Calling load_response"
    response_function = load_response(data, energy_range)

    @mpiinfo "Generating emissions model"
    emission_model = prepare_model_mekal(nHcol, energy_range, redshift, use_interpolation=use_interpolation)

    @mpiinfo "Testing emissions model"
    em_direct = prepare_model_mekal(nHcol, energy_range, redshift, use_interpolation=false)

    model = Model_NFW_GNFW(
        5e14u"Msun",
        0.13,
        1.0510, # Using universal values from Arnaud 2010
        5.4905,
        0.3081,
        1.177,
        redshift,
        [64, 64],
        0.492u"arcsecondᵃ",
        emission_model,
        100e3u"s",
        response_function,
        (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        center_radius
    )
    model_direct = Model_NFW_GNFW(
        5e14u"Msun",
        0.13,
        1.0510, # Using universal values from Arnaud 2010
        5.4905,
        0.3081,
        1.177,
        redshift,
        [64, 64],
        0.492u"arcsecondᵃ",
        em_direct,
        100e3u"s",
        response_function,
        (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        center_radius
    )
    replace!(model, NaN => 0)
    replace!(model_direct, NaN => 0)

    err = sum(abs2, model - model_direct)
    @mpiinfo "Error in emission model is" err

    sample(
        obs,
        bg,
        response_function,
        prior_set,
        observation.second,
        observed_background.second,
        redshift;
        emission_model=emission_model,
        pixel_edge_angle=pixel_edge_angle,
        center_radius=center_radius
    )
end
