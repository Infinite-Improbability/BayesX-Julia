include("likelihood.jl")

struct Results{T<:Real,A<:AbstractArray{T},V<:AbstractVector{Prior},R<:Real}
    observation::A
    background::A
    response_function::Matrix
    cluster_model::Function
    priors::V
    nHcol::SurfaceDensity
    redshift::Real
    x::NTuple{2,<:Real}
    y::NTuple{2,<:Real}
    bin_size::Real
    use_interpolation::Bool
    centre_radius::R
    mask::Union{Matrix{Bool},Nothing}
end