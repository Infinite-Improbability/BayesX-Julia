struct Results
    observation::Array{<:T}
    background::Array{<:T}
    response_function::Matrix
    cluster_model::Function
    priors <: AbstractVector{<:Prior}
    nHcol::SurfaceDensity
    redshift::Real,
    x::NTuple{2,<:Real},
    y::NTuple{2,<:Real}
    bin_size::Real,
    use_interpolation::Bool,
    centre_radius <: Real,
    mask::Union{Matrix{Bool},Nothing}
end