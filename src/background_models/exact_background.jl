using PoissonRandom
export Exact_Background

"""
    Exact_Background(background_counts::Array{Int64}; scale::Float64, kwargs...)

    Takes an exact background grid (channel, x, y) or (channel) and scales it by a factor. 
"""
function Exact_Background(background_counts::AbstractArray, scale::Float64, kwargs...)
    @assert scale >= 0
    return background_counts .* scale
end