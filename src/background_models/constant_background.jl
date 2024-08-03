using PoissonRandom
export Uniform_Background

"""
    Uniform_Background(background_rate:Number, num_channels:Integer; kwargs...)

Creates a poisson sampled uniform energy spectrum for the background. Intended for testing purposes.
"""
function Uniform_Background(background_rate::Number, num_channels::Integer; kwargs...)
    @assert background_rate >= 0

    background = Array{Float64}(undef, num_channels)
    for i in eachindex(background)
        background[i] = background_rate
    end
    return background
end