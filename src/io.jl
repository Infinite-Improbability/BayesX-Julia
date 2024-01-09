using Unitful, DimensionfulAngles
using ProgressMeter

"""
Abstract type to wrap various specific formats for source data. By doing so we can make use of
    multiple dispatch to automatically select the appropriate loading functions for data.
"""
abstract type Dataset end
# Include exposure times in dataset
# And NHcol, bg count rate?

include("io/fits.jl")
include("io/plaintext.jl")

"""
    load_data(data)

Load events data from a given dataset, returning pairs of events and exposure times.

The first pair includes all observed events, the second is the background observation.
"""
function load_data(data::Dataset)::NTuple{2,Pair}
    error("load_data is not implemented for these arguments", data)
end

"""
    load_response(data, min_energy, max_energy)

Load the RMF and ARF for an observation, trimmed for the desired energy range
"""
function load_response(data::Dataset, min_energy::Unitful.Energy, max_energy::Unitful.Energy)
    error("load_response is not implemented for these arguments\n", "data = $data\n", "energy_range = $min_energy,$max_energy")
end


"""
    load_mask(data)

The matrix is true for masked pixels.
"""
function load_mask(data::Dataset)::Matrix{Bool}
    error("load_data is not implemented for these arguments", data)
end

"""
    bin_events(::Dataset, events, ...)

Placeholder function, returns events matrix without modification"""
function bin_events(::Dataset, events, args...)
    @mpiwarn "bin_events is not implemented for these arguments."
    return events
end

"""
    Ellipse(x, y, r1, r2, θ)

An ellipse centered at `(x, y)` with radii `r1` and `r2` and rotation `θ`.

`r1` and `r2` are the width and height respectively with zero rotation and
`θ` is defined as rotation counterclockwise from the x-axis in radians.
"""
struct Ellipse{T<:AbstractFloat}
    x::T
    y::T
    r1::T
    r2::T
    θ::T
end

"""
    test_point(e::Ellipse, x, y; atol=0)

Tests whether a point is within an ellipse with tolerance `atol`.

The equation for an ellipse in Cartesian coordinates is of the form `f(x,y) = 1`. We consider
any point `(x0, y0)` such that `f(x0, y0) <= 1 + atol` as being within the ellipse.
"""
function test_point(e::Ellipse, x, y; atol=0)::Bool
    t1 = (cos(e.θ) * (x - e.x) + sin(e.θ) * (y - e.y))^2
    t2 = (sin(e.θ) * (x - e.x) - cos(e.θ) * (y - e.y))^2

    return (t1 / (e.r1^2) + t2 / (e.r2^2)) <= (1 + atol)
end
function test_point(e::Ellipse, xy::AbstractVector; atol=0)::Bool
    return test_point(e, xy...; atol=atol)
end

"""
    load_mask(path, x_edges, y_edges)

Takes a mask file as used by DS9 and x and y edges to bin on. This should match
the edges used to bin the data.
"""
function load_mask(path::AbstractString, x_edges, y_edges)::Matrix{Bool}
    ellipses = Set{Ellipse}()
    # Convert lines of ellipse definitions in input file to a list of Ellipse objects.
    open(path) do f
        for line in readlines(f)
            if !occursin("ellipse", line)
                @mpidebug "Skipping mask line" line
                continue
            end
            @mpidebug "Processing mask line" line
            # Currently assuming pixel coordinates. In theory this isn't guaranteed.
            # format is ellipse(x, y, r1, 2, angle)
            # first strip ellipse( to get x, y, r1, 2, angle)
            line = split(line, '(')[end]
            # then strip ) to get x, y, r1, 2, angle
            line = rstrip(line, ')')
            # then split on commas and cast to floats
            args = Tuple(parse(Float64, el) for el in split(line, ","))
            push!(ellipses, Ellipse(args...))
        end
    end

    @mpidebug "Loaded mask ellipses" length(ellipses)

    centres_x = (x_edges[1:end-1] + x_edges[2:end]) / 2
    centres_y = (y_edges[1:end-1] + y_edges[2:end]) / 2

    mask = zeros(Bool, length(centres_x), length(centres_y))
    pixels_masked = 0

    for j in eachindex(centres_x)
        for i in eachindex(centres_y)
            mask[i, j] = any(test_point.(ellipses, x_edges[i], y_edges[j]))
            if mask[i, j]
                pixels_masked += 1
            end
        end
    end

    @mpiinfo "Mask loaded" pixels_masked

    return mask
end