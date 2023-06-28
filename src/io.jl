using FITSIO

# Import our patch
include("fitsio_fix.jl")

struct BayesXData{N<:Int,T<:AbstractFloat}
    observation::Matrix{N}
    background::Matrix{N}
    arf::Vector{T}
    rmf::Matrix{T}
end

function load_data(events::AbstractString, background::AbstractString)
end

function safe_read_key(hdu::HDU, key::String)
    try
        return read_key(hdu, key)
    catch
        @debug "Exception when looking for events file. Probably a HDU without extname."
        return (nothing, nothing)
    end

end

function load_fits_file(path::AbstractString)#::Matrix{Int64}
    f = FITS(path, "r")
    event_hdus::Vector{TableHDU} = [h for h in f if safe_read_key(h, "extname")[1] == "EVENTS"]
    if length(event_hdus) > 1
        @warn "$(length(event_hdus)) HDUs with events extension found. Using the first."
    end
    h = event_hdus[1]
    @info "Selected HDU with HDUNAME '$(safe_read_key(h, "HDUNAME")[1])'"
end