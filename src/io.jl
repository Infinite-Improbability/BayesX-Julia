using FITSIO

# Import our patch
include("fitsio_fix.jl")

abstract type BayesXDataset end
# Include exposure times in dataset
# And NHcol, bg count rate?

struct FITSData{S<:AbstractString} <: BayesXDataset
    observation::S
    background::S
    arf::S
    rmf::S
end

function load_data(data::BayesXDataset)
end

function load_response(data::FITSData)
end

function safe_read_key(hdu::HDU, key::String, msg::AbstractString)
    try
        return read_key(hdu, key)
    catch
        @debug msg
        return (nothing, nothing)
    end

end

function load_events_from_fits(path::AbstractString)
    f = FITS(path, "r")
    event_hdus::Vector{TableHDU} = [
        h for h in f
        if safe_read_key(h, "extname", "Exception when looking for events HDU. Probably a HDU without extname.")[1] == "EVENTS"
    ]
    if length(event_hdus) > 1
        @warn "$(length(event_hdus)) HDUs with events extension found. Using the first."
    end
    h = event_hdus[1]
    @info "Selected HDU with HDUNAME '$(read_key(h, "HDUNAME")[1])'"

    return [read(h, "x") read(h, "y") read(h, "pi")]
end

function load_data(data::FITSData)
    obs = load_events_from_fits(data.observation)
    bg = load_events_from_fits(data.background)
    return (obs, bg)
end

function load_response(data::FITSData)
    f = FITS(data.rmf)
    rmf_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname", "Exception when looking for matrix HDU. Probably a HDU without extname.")[1] == "MATRIX"]
    if length(matrix_hdus) > 1
        @warn "$(length(matrix_hdus)) HDUs with matrix extension found. Using the first."
    end
    r = rmf_hdus[1]
    @info "Selected HDU with HDUNAME '$(read_key(r, "HDUNAME")[1])'"

    first_channel = read(r, "F_CHAN")
    last_channel = first_channel .+ read(r, "N_CHAN") - 1
    channels_for_bin = read(r, "MATRIX")

    rmf = zeros(Float64, (maximum(last_channel), read_key(r, "NAXIS2")))

    for i in axes(rmf, 2)
        rmf[first_channel[i]:last_channel[i], i] .= channels_for_bin[i]
    end

    f = FITS(data.rmf)
    arf_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname", "Exception when looking for specrsp HDU. Probably a HDU without extname.")[1] == "MATRIX"]
    if length(matrix_hdus) > 1
        @warn "$(length(matrix_hdus)) HDUs with matrix extension found. Using the first."
    end
    a = arf_hdus[1]
    @info "Selected HDU with HDUNAME '$(read_key(r, "HDUNAME")[1])'"
    arf = read(a, "SPECRESP")

    # turn rmf into response matrix
    for i in axes(rmf, 1)
        rmf[i, :] .*= arf
    end

    return rmf
end

function bin_events(events, energy_range, x_edges, y_edges)::Array{Int64}
    events = events[minimum(x_edges).<events[:, 1].<maximum(x_edges), :]
    events = events[minimum(y_edges).<events[:, 2].<maximum(y_edges), :]
    events = events[energy_range[1].<events[:, 3].<energy_range[2], :]

    binned = zeros(Int64, (length(x_edges) - 1, length(y_edges) - 1, trunc(Int64, maximum(events[:, 3]))))

    for r in eachrow(events)
        i = searchsortedlast(x_edges, r[1])
        j = searchsortedlast(y_edges, r[2])
        binned[i, j, trunc(Int64, r[3])] += 1
    end

    return binned
end