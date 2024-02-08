using Unitful, UnitfulAstro
using ArgCheck

export Model_Piecewise

function validate_interpolation_args(x::Vector{T}, y::Vector{U})::Bool where {T<:Number,U<:Number}
    return length(x) == length(y) && all(diff(x) .> zero(T))
end

"""
    LinearInterpolator(x, y)

Construct a linear interpolator for the given data.

# Boundary behaviour

If `x` is below the minimum value of `x` in the data it is rounded up the the minimum. This
results in a flat line below the minimum value of `x`. In the cluster context this means the
core properties may be constant below a certain radius. To avoid this simply fix the smallest value
of `x`` at zero.

If `x` is above the maximum value of `x` in the data the interpolator returns zero. This is
because the density and temperature are assumed to be zero outside the cluster. Fix the density
and temperature at your largest radius to zero if you wish to avoid a discontinuity here.
"""
struct LinearInterpolator{T<:Number,U<:Number}
    x::Vector{T}
    y::Vector{U}
    min_x::T
    max_x::T

    LinearInterpolator(x::Vector{T}, y::Vector{U}) where {T<:Number,U<:Number} =
        validate_interpolation_args(x, y) ? new{T,U}(x, y, extrema(x)...) : throw(ArgumentError("Arguments must have same length."))
end
function (li::LinearInterpolator)(x::Number)
    if x < li.min_x
        x = li.min_x
    elseif x > li.max_x
        return zero(eltype(li.y))
    elseif x == li.max_x
        return li.y[end] # the lower bound would be fine but the upper bound would trigger index errors
    end

    i = searchsortedlast(li.x, x)

    return li.y[i] + (li.y[i+1] - li.y[i]) / (li.x[i+1] - li.x[i]) * (x - li.x[i])
end

"""
    Model_Piecewise(radius, density, temperature, radius, density, temperature, ...)

Construct a piecewise model for the given data.

The model is constructed by linearly interpolating the given data. See [`LinearInterpolator`](@ref)
for details on the interpolation behaviour.

Density is assumed to be in units of `g/cm^3`, temperature in units of `keV`, and radius in units of `kpc`.
"""
function Model_Piecewise(args::Vararg{Real}; kwargs...)
    @argcheck length(args) % 3 == 0 "Model_Piecewise args must be a multiple of three, in a sequence `r_i`, `ρ_i`, `T_i`"

    radii = collect(args[1:3:end]) * 1u"kpc"
    densities = collect(args[2:3:end]) * 1u"g/cm^3"
    temperatures = collect(args[3:3:end]) * 1u"keV"

    @argcheck all(radii .>= 0u"kpc") "Model_Piecewise radii must be non-negative"
    @argcheck all(densities .>= 0u"g/cm^3") "Model_Piecewise densities must be non-negative"
    @argcheck all(temperatures .>= 0u"keV") "Model_Piecewise temperatures must non-negative"

    sorted = sortperm(radii)
    radii = radii[sorted]
    densities = densities[sorted]
    temperatures = temperatures[sorted]

    density_interpolator = LinearInterpolator(radii, densities)
    temperature_interpolator = LinearInterpolator(radii, temperatures)

    return temperature_interpolator, density_interpolator
end