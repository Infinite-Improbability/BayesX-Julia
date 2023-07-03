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

function load_data(data::BayesXData, energy_range)::Tuple{Matrix}
end

function safe_read_key(hdu::HDU, key::String)
    try
        return read_key(hdu, key)
    catch
        @debug "Exception when looking for events file. Probably a HDU without extname."
        return (nothing, nothing)
    end

end

function load_events_from_fits(path::AbstractString)
    f = FITS(path, "r")
    event_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname")[1] == "EVENTS"]
    if length(event_hdus) > 1
        @warn "$(length(event_hdus)) HDUs with events extension found. Using the first."
    end
    h = event_hdus[1]
    @info "Selected HDU with HDUNAME '$(safe_read_key(h, "HDUNAME")[1])'"

    return [h["X"] h["y"] h["PI"]]
end

function load_data(data::FITSData, energy_range)
    obs = load_events_from_fits(data.observation)
end

