using Unitful

# This should be replaced with proper background modelling
function prepare_background(
    observed::AbstractArray{<:Integer,3},
    observed_background::AbstractArray{<:Integer,3},
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
)::NTuple{2,Vector{Float64}}
    # implicitly includes average effective area and pixel edge angle
    bg_count_rate = [mean(@view observed_background[i, :, :]) for i in axes(observed_background, 1)] ./ bg_exposure_time
    zero_channels = count(i -> i == 0u"s^-1", bg_count_rate)
    total_channels = size(observed_background, 1)
    fallback_count_rate = 1.0e-7 / bg_exposure_time
    if zero_channels > 0
        @mpiwarn "Some energy channels have zero background counts. Replacing with fallback value" zero_channels total_channels fallback_count_rate
        replace!(bg_count_rate, 0u"s^-1" => fallback_count_rate)
    end

    @assert all(i -> i > 0u"s^-1", bg_count_rate)
    @mpidebug "Background rate estimated" bg_count_rate

    # vector of background as a function of energy
    # counts per pixel per exposure time
    predicted_obs_bg = bg_count_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = bg_count_rate * bg_exposure_time # Used for log likelihood

    @assert all(isfinite, predicted_obs_bg)
    @assert all(isfinite, predicted_bg_bg)
    @assert all(i -> i > 0, predicted_obs_bg)
    @assert all(i -> i > 0, predicted_bg_bg)
    @assert length(predicted_bg_bg) == size(observed, 1)

    return predicted_obs_bg, predicted_bg_bg
end