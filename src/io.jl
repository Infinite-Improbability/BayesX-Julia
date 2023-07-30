using FITSIO
using Unitful, DimensionfulAngles
using ProgressMeter
using DelimitedFiles

export FITSData, PlaintextData

# Import our patch
include("fitsio_fix.jl")

"""
Abstract type to wrap various specific formats for source data. By doing so we can make use of
    multiple dispatch to automatically select the appropriate loading functions for data.
"""
abstract type Dataset end
# Include exposure times in dataset
# And NHcol, bg count rate?

# TODO: Be more specific with types
"""
Store a collection of observation data from fits files.
"""
struct FITSData{S<:AbstractString,T<:DimensionfulAngles.Angle} <: Dataset
    observation::S
    background::S
    arf::S
    rmf::S
    pixel_edge_angle::T
end

"""
    load_data(data)

Load events data from a given dataset, returning pairs of events and exposure times.

The first pair includes all observed events, the second is the background observation.
"""
function load_data(data::Dataset)::NTuple{2,Pair}
    error("load_data is not implemented for these arguments", data)
end

"""
    load_response(data, energy_range)

Load the RMF and ARF for an observation, trimmed for the desired energy range
"""
function load_response(data::Dataset, energy_range)::Matrix{Float64}
    error("load_response is not implemented for these arguments", data, energy_range)
end

function safe_read_key(hdu::HDU, key::String, msg::AbstractString)
    try
        return read_key(hdu, key)
    catch
        @mpidebug msg
        return (nothing, nothing)
    end

end

"""
    load_events_from_fits(path)

Loads events from a single fits file for further processing.
"""
function load_events_from_fits(path::AbstractString)::Pair{Matrix,Unitful.Time{Float64}}
    f = FITS(path, "r")
    event_hdus::Vector{TableHDU} = [
        h for h in f
        if safe_read_key(h, "extname", "Exception when looking for events HDU. Probably a HDU without extname.")[1] == "EVENTS"
    ]
    if length(event_hdus) > 1
        @warn "$(length(event_hdus)) HDUs with events extension found. Using the first."
    end
    h = event_hdus[1]
    @mpidebug "Selected HDU with EVENTS extension and HDUNAME '$(safe_read_key(h, "HDUNAME", "HDU has no name")[1])' from $path for events."

    live_time = read_key(h, "LIVETIME")[1] # We don't want safe_read_key because we want an exception if this fails.

    return [read(h, "x") read(h, "y") read(h, "pi") read(h, "energy") * 1u"eV"] => live_time * 1u"s" # TODO: manual units bad
end

"""
    load_data(data::FITSData)

Load data from FITS events files for the observation and background.

Exposure times are taken from the LIVETIME value in the header.
"""
function load_data(data::FITSData)::NTuple{2,Pair}
    obs = load_events_from_fits(data.observation)
    bg = load_events_from_fits(data.background)
    return (obs, bg)
end

function load_response(data::FITSData, energy_range)::Matrix{Unitful.Area{Float64}}
    @mpidebug "Loading response matrices"

    f_rmf = FITS(data.rmf)
    rmf_hdus::Vector{TableHDU} = [
        h for h in f_rmf if safe_read_key(h, "extname", "Exception when looking for matrix HDU. Probably a HDU without extname.")[1] == "MATRIX"
    ]
    if length(rmf_hdus) > 1
        @warn "$(length(rmf_hdus)) HDUs with matrix extension found. Using the first."
    end
    r = rmf_hdus[1]
    @mpidebug "Selected HDU with MATRIX extension and HDUNAME '$(safe_read_key(r, "HDUNAME", "HDU has no name")[1])' from $(data.rmf) for RMF"

    first_channel = [i[1] for i in read(r, "F_CHAN")]
    last_channel = first_channel .+ [i[1] for i in read(r, "N_CHAN")] .- 1
    channels_for_bin = read(r, "MATRIX")

    rmf = zeros(Float64, (maximum(last_channel), read_key(r, "NAXIS2")[1]))

    for i in axes(rmf, 2)
        rmf[first_channel[i][1]:last_channel[i][1], i] .= channels_for_bin[i]
    end

    f_arf = FITS(data.arf)
    arf_hdus::Vector{TableHDU} = [
        h for h in f_arf if safe_read_key(h, "extname", "Exception when looking for specrsp HDU. Probably a HDU without extname.")[1] == "SPECRESP"
    ]
    if length(arf_hdus) > 1
        @warn "$(length(arf_hdus)) HDUs with matrix extension found. Using the first."
    end
    a = arf_hdus[1]
    @mpidebug "Selected HDU with SPECRESP extension and HDUNAME '$(safe_read_key(r, "HDUNAME", "HDU has no name")[1])' from $(data.arf) for ARF"
    arf = read(a, "SPECRESP") # TODO: Automate unit selection

    # turn rmf into response matrix
    # now our variable name is inaccurate but we don't have to allocate new memory
    for i in axes(rmf, 1)
        rmf[i, :] .*= arf
    end

    # Get first bin where maximum energy >= min of range
    # ∴ Subsequent bins must have minimum energy > min of range
    min_bin = searchsortedfirst(read(r, "ENERG_HI") * 1u"keV", minimum(energy_range))
    # min_channel = searchsortedfirst(read(f_rmf[3], "E_MAX") * 1u"keV", minimum(energy_range))
    min_channel = 1 # we don't trim events channels based on minimum (yet)
    # Get last bin where minimum energy >= max of range
    # ∴ This and subsequent bins must have minimum energy > max of range
    max_bin = searchsortedfirst(read(r, "ENERG_LO") * 1u"keV", maximum(energy_range)) - 1
    max_channel = searchsortedfirst(read(f_rmf[3], "E_MIN") * 1u"keV", maximum(energy_range)) - 1

    @mpidebug "Trimming response matrix with arrangement (PI,E) to range" min_channel max_channel min_bin max_bin

    return rmf[min_channel:max_channel, min_bin:max_bin] * 1u"cm^2" # Hack to add arf units
end

"""
    bin_events(::FITSData, events, energy_range, x_edges, y_edges)

Take a table of events in the format (x y channel energy), trim it by energy and bin it spatially.
Returns an array of counts per bin with dimensions (channel, x, y).
"""
function bin_events(_::FITSData, events, energy_range, x_edges, y_edges)::Array{Int64}
    @mpidebug "Binning events"

    events = events[minimum(x_edges).<events[:, 1].<maximum(x_edges), :]
    events = events[minimum(y_edges).<events[:, 2].<maximum(y_edges), :]
    events = events[minimum(energy_range).<events[:, 4].<maximum(energy_range), :]

    # Maximum doesn't like the Unitful array
    channels = ustrip.(events[:, 3])

    binned = zeros(Int64, (trunc(Int64, maximum(channels)), length(x_edges) - 1, length(y_edges) - 1))

    prog = ProgressUnknown("Events read:")
    for r in eachrow(events)
        i = searchsortedlast(x_edges, r[1])
        j = searchsortedlast(y_edges, r[2])
        binned[trunc(Int64, r[3]), i, j] += 1
        ProgressMeter.next!(prog)
    end
    ProgressMeter.finish!(prog)

    return binned
end

"""Store preprocessed data from plaintext files with whitespace delimiters."""
struct PlaintextData{S<:AbstractString,T<:Unitful.Time,U<:DimensionfulAngles.Angle,V<:Unitful.Area} <: Dataset
    observation::S
    background::S
    arf::S
    rmf::S
    shape::NTuple{3,Int}
    rmf_shape::NTuple{2,Int}
    obs_exposure_time::T
    bg_exposure_time::T
    pixel_edge_angle::U
    rmf_unit::V
end

"""
    load_data(data::PlaintextData)

Load data from plaintext files as generated by existing Bayes-X binning scripts.

Data is reshaped to the given size so should be provided in column-major order
or in a matrix of correct size.
"""
function load_data(data::PlaintextData)::NTuple{2,Pair}
    obs_data = reshape(readdlm(data.observation, Int), data.shape)
    bg_data = reshape(readdlm(data.background, Int), data.shape)
    return (obs_data => data.obs_exposure_time, bg_data => data.bg_exposure_time)
end

"""
    load_response(data::PlaintextData)

Load response matrix from plaintext files as generated by existing Bayes-X binning scripts.

Data is reshaped to the given size so should be provided in column-major order
or as a matrix of correct size.
"""
function load_response(data::PlaintextData)::NTuple{2,Pair}
    rmf = reshape(readdlm(data.rmf), data.rmf_shape)
    arf = reshape(readdlm(data.arf), data.rmf_shape[2])

    for i in axes(rmf, 1)
        rmf[i, :] .*= arf
    end

    return rmf * data.unit
end

"""
    bin_events(::PlaintextData, events, ...)

Placeholder function, returns events matrix without modification"""
function bin_events(_::PlaintextData, events, args...)
    @mpiwarn "PlaintextData does not support binning. Data will be used unmodified."
    return events
end