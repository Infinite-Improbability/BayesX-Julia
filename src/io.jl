using FITSIO

# Import our patch
include("fitsio_fix.jl")

abstract type BayesXDataset end
# Include exposure times in dataset
# And NHcol, bg count rate?

# TODO: Be more specific with types
struct FITSData{S<:AbstractString} <: BayesXDataset
    observation::S
    background::S
    arf::S
    rmf::S
    exposure_time # TODO: Decouple source and bg exposure times
    pixel_edge_angle
end

function load_data(data::BayesXDataset)
end

function load_response(data::BayesXDataset, energy_range)
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
    @info "Selected HDU with HDUNAME '$(safe_read_key(h, "HDUNAME", "HDU has no name")[1])'"

    return [read(h, "x") read(h, "y") read(h, "pi") read(h, "energy") * 1u"eV"] # TODO: manual units bad
end

function load_data(data::FITSData)
    obs = load_events_from_fits(data.observation)
    bg = load_events_from_fits(data.background)
    return (obs, bg)
end

function load_response(data::FITSData, energy_range)
    f = FITS(data.rmf)
    rmf_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname", "Exception when looking for matrix HDU. Probably a HDU without extname.")[1] == "MATRIX"]
    if length(rmf_hdus) > 1
        @warn "$(length(rmf_hdus)) HDUs with matrix extension found. Using the first."
    end
    r = rmf_hdus[1]
    @info "Selected HDU with HDUNAME '$(safe_read_key(r, "HDUNAME", "HDU has no name")[1])'"

    first_channel = [i[1] for i in read(r, "F_CHAN")]
    last_channel = first_channel .+ [i[1] for i in read(r, "N_CHAN")] .- 1
    channels_for_bin = read(r, "MATRIX")

    rmf = zeros(Float64, (maximum(last_channel), read_key(r, "NAXIS2")[1]))

    for i in axes(rmf, 2)
        rmf[first_channel[i][1]:last_channel[i][1], i] .= channels_for_bin[i]
    end

    f = FITS(data.arf)
    arf_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname", "Exception when looking for specrsp HDU. Probably a HDU without extname.")[1] == "SPECRESP"]
    if length(arf_hdus) > 1
        @warn "$(length(arf_hdus)) HDUs with matrix extension found. Using the first."
    end
    a = arf_hdus[1]
    @info "Selected HDU with HDUNAME '$(safe_read_key(r, "HDUNAME", "HDU has no name")[1])'"
    arf = read(a, "SPECRESP") # TODO: Automate unit selection

    # turn rmf into response matrix
    for i in axes(rmf, 1)
        rmf[i, :] .*= arf
    end

    # Get first bin where maximum energy >= min of range
    # ∴ Subsequent bins must have minimum energy > min of range
    min_bin = searchsortedfirst(read(r, "ENERG_HI") * 1u"keV", energy_range[1])
    min_channel = first_channel[min_bin]
    # Get last bin where minimum energy >= max of range
    # ∴ This and subsequent bins must have minimum energy > max of range
    max_bin = searchsortedfirst(read(r, "ENERG_LO") * 1u"keV", energy_range[2]) - 1
    max_channel = last_channel[max_bin] # assumes max bin increases montonically

    return rmf[:, min_bin:max_bin] * 1u"cm^2" # Hack to add arf units
end

function bin_events(events, energy_range, x_edges, y_edges)::Array{Int64}
    events = events[minimum(x_edges).<events[:, 1].<maximum(x_edges), :]
    events = events[minimum(y_edges).<events[:, 2].<maximum(y_edges), :]
    events = events[energy_range[1].<events[:, 4].<energy_range[2], :]

    # display(events)

    # Maximum doesn't like the Unitful array
    channels = ustrip.(events[:, 3])

    binned = zeros(Int64, (length(x_edges) - 1, length(y_edges) - 1, trunc(Int64, maximum(channels))))

    for r in eachrow(events)
        i = searchsortedlast(x_edges, r[1])
        j = searchsortedlast(y_edges, r[2])
        binned[i, j, trunc(Int64, r[3])] += 1
    end

    return transpose(binned)
end