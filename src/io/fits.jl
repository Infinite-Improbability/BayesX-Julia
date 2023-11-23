# Import our patch
using FITSIO

export FITSData

include("fitsio_fix.jl")

# TODO: Be more specific with types

"""
    FITSData(observation::AbstractString, background::AbstractString, arf::AbstractString, rmf::AbstractString, pixel_edge_angle::DimensionfulAngles.Angle)

Store a collection of observation data from fits files.

Takes paths to the files as inputs.
"""
struct FITSData{S<:AbstractString,T<:DimensionfulAngles.Angle} <: Dataset
    observation::S
    background::S
    arf::S
    rmf::S
    pixel_edge_angle::T
end
function FITSData(obs::S, bg::S, arf::S, rmf::S, pea::Unitful.DimensionlessQuantity) where {S<:AbstractString}
    @warn "You should use dimensionful angles"
    FITSData(obs, bg, arf, rmf, ustrip(u"rad", pea) * 1u"radᵃ")
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
function load_events_from_fits(path::AbstractString)::Pair{<:Matrix,<:Unitful.Time{Float64}}
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

    live_time = safe_read_key(h, "LIVETIME", "LIVETIME key absent, trying EXPOSURE")[1]
    if isnothing(live_time)
        live_time = read_key(h, "EXPOSURE")[1] # We don't want safe_read_key because we want an exception if this fails.
    end

    return [read(h, "x") read(h, "y") read(h, "pi")] => live_time * 1u"s" # TODO: manual units bad
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

"""
    load_response(data::FITSData, min_energy::Unitful.Energy, max_energy::Unitful.Energy)

Load response function for specified energy range, adjusting it to align with channel edges.

Returns response matrix, energy bins and included channels.
"""
function load_response(data::FITSData, min_energy::Unitful.Energy, max_energy::Unitful.Energy)::Tuple{Matrix{<:Unitful.Area{Float64}},Vector{<:Unitful.Energy},NTuple{2,Int64}}
    @mpidebug "Loading response matrices"

    # TODO: Automate units
    @mpidebug "Assuming keV and cm^2 for RMF and ARF units"

    #### RMF ####
    # The RMF is stored in the format
    # ENERG_LO ENERG_HI N_GRP F_CHAN N_CHAN MATRIX

    # Open RMF file
    f_rmf = FITS(data.rmf)
    rmf_hdus::Vector{TableHDU} = [
        h for h in f_rmf if safe_read_key(h, "extname", "Exception when looking for matrix HDU. Probably a HDU without extname.")[1] == "MATRIX"
    ]
    if length(rmf_hdus) > 1
        @warn "$(length(rmf_hdus)) HDUs with matrix extension found. Using the first."
    end
    rmf = rmf_hdus[1]
    @mpidebug "Selected HDU with MATRIX extension and HDUNAME '$(safe_read_key(rmf, "HDUNAME", "HDU has no name")[1])' from $(data.rmf) for RMF"

    # ENERG_LO and ENERGY_HI are the bounding energies of the bins.
    # We verify these bins are continuous and not overlapping.
    rmf_energy_bin_minimums = read(rmf, "ENERG_LO") * 1u"keV"
    rmf_energy_bin_maximums = read(rmf, "ENERG_HI") * 1u"keV"
    @assert all(rmf_energy_bin_minimums[2:end] .== rmf_energy_bin_maximums[1:end-1]) "The energy bins of the RMF are not continuous."
    for i in axes(rmf_energy_bin_minimums, 1)[1:end-1]
        @assert rmf_energy_bin_minimums[i+1] > rmf_energy_bin_minimums[i] "The RMF energy bins are not monotonically increaing."
    end

    # Get channel ranges for each energy bin
    # F_CHAN is the first channel for a bin
    # N_CHAN is the number of channels for a bin
    # For mysterious reasons both are given as vectors of length one
    # We'll verify that in case anyone ever passes in a RMF that behaves differently
    @assert all(length.(read(rmf, "F_CHAN")) .== 1) "Multiple first channel (F_CHAN) values found for single energy bin."
    @assert all(length.(read(rmf, "N_CHAN")) .== 1) "Multiple number of channels (N_CHAN) values found for single energy bin."
    # Then we read in the first channel for each energy bin
    first_channel_of_bin = [i[1] for i in read(rmf, "F_CHAN")]
    # And the number of channels
    number_of_channels_in_bin = [i[1] for i in read(rmf, "N_CHAN")]
    # And calculate the last channel as first_channel + number_of_channels - 1
    # The -1 accounts for the first channel channel contributing to the total number of channels
    last_channel_of_bin = first_channel_of_bin .+ number_of_channels_in_bin .- 1

    # Then we load in MATRIX
    # For each bin this a vector giving the probability of a photon going to each channel
    # So we want the vector to be the same length as the number of channels in the bin
    energy_to_channel_mapping = read(rmf, "MATRIX")
    @assert length.(energy_to_channel_mapping) == number_of_channels_in_bin "Discrepency between stated number of channels in bin (N_CHAN) and length of energy to channel mapping vector (MATRIX)"

    # Preallocate output matrix
    # Response matrix has shape (PI, E) where PI is channel and E is energy bin
    # NAXIS gives the number of rows (energy bins) in the table
    # We'll verify this matches our other data
    @assert read_key(rmf, "NAXIS2")[1] == length(first_channel_of_bin) "Claimed number of rows in table does not match length of data"
    response_matrix = zeros(Float64, (maximum(last_channel_of_bin), read_key(rmf, "NAXIS2")[1]))

    # tidy up
    close(f_rmf)

    # Load values into matrix
    # We iterate over energy bins, which are our columns
    @views for bin in axes(response_matrix, 2)
        # for each bin `i` we select the appropriate channel range from the matrix with first_channel_of_bin[i]:last_channel_of_bin[i]
        # and assign the energy mapping to it
        response_matrix[first_channel_of_bin[bin]:last_channel_of_bin[bin], bin] = energy_to_channel_mapping[bin]
    end

    # The RMF gives the probability of a photon going into each channel if it falls into a given energy bin
    # Thus the sum of each bin should be <= 1
    # We'll assume every photon is assigned to a channel
    # Then we can assert the sum ≈ 1
    # In practice it appears the sum can sometimes be slightly above one
    # I'm attributing this to floating point errors
    col_sum = sum(response_matrix, dims=1)
    @assert all(@. isapprox(col_sum, 1.0, rtol=0.001) || col_sum < 1) "Invalid RMF column sum"

    # Now we have the RMF in matrix form RMF(PI, E)

    #### ARF ####
    # The ARF is stored in the format
    # ENERG_LO ENERG_HI SPECRESP

    # Open ARF file
    f_arf = FITS(data.arf)
    arf_hdus::Vector{TableHDU} = [
        h for h in f_arf if safe_read_key(h, "extname", "Exception when looking for specrsp HDU. Probably a HDU without extname.")[1] == "SPECRESP"
    ]
    if length(arf_hdus) > 1
        @warn "$(length(arf_hdus)) HDUs with matrix extension found. Using the first."
    end
    arf = arf_hdus[1]
    @mpidebug "Selected HDU with SPECRESP extension and HDUNAME '$(safe_read_key(rmf, "HDUNAME", "HDU has no name")[1])' from $(data.arf) for ARF"

    # Verify the energy bins match the RMF
    @assert rmf_energy_bin_minimums == (read(arf, "ENERG_LO") * 1u"keV") "Lower edges of ARF energy bins do not match lower edge for RMF energy bins."
    @assert rmf_energy_bin_maximums == (read(arf, "ENERG_HI") * 1u"keV") "Upper edges of ARF energy bins do not match Upper edge for RMF energy bins."

    # Load effective area per bin
    effective_area_per_energy_bin = read(arf, "SPECRESP")
    @assert length(effective_area_per_energy_bin) == size(response_matrix, 2) "ARF has a different number of energy bins to RMF"

    # tidy up
    close(f_arf)

    # Turn rmf into response matrix by multiplying with arf
    # Each bin column is multiplied by the appropriate effective area.
    for bin in axes(response_matrix, 2)
        response_matrix[:, bin] .*= effective_area_per_energy_bin[bin]
    end

    # Get first bin where lower edge >= user specified `min_energy``
    min_bin = searchsortedfirst(rmf_energy_bin_minimums, min_energy)
    min_channel = first_channel_of_bin[min_bin]

    # Get last bin where upper edge <= user specificed `max_energy``
    max_bin = searchsortedlast(rmf_energy_bin_maximums, max_energy)
    max_channel = last_channel_of_bin[max_bin]

    @mpidebug "Trimming response matrix with arrangement (PI,E) to range" min_channel max_channel min_bin max_bin
    trimmed_response_matrix = response_matrix[min_channel:max_channel, min_bin:max_bin] * 1u"cm^2" # we only apply units here for type reasons

    # Now we want an energy range
    # Convert bounding bin indices to energy values
    new_min_energy = rmf_energy_bin_minimums[min_bin]
    new_max_energy = rmf_energy_bin_maximums[max_bin]

    # Trim bin edges to the selected bin range
    # This could be merged with the next step but this helps readability
    bin_mins = rmf_energy_bin_minimums[min_bin:max_bin]
    bin_maxes = rmf_energy_bin_maximums[min_bin:max_bin]

    # Load in all the bin mins and add the last max to end the upper bin
    # Now we have a vector of energy bin edges
    energy_bins = [bin_mins; bin_maxes[end]] # note the semicolon to concatenate
    n_energy_bins = length(energy_bins)

    @mpiinfo "Energy range adjusted to align with response matrix" new_min_energy new_max_energy n_energy_bins

    return trimmed_response_matrix, energy_bins, (min_channel, max_channel)
end

"""
    bin_events(::FITSData, events, energy_range, x_edges, y_edges)

Take a table of events in the format (x y channel), trim it by energy and bin it spatially.
Returns an array of counts per bin with dimensions (channel, x, y).
"""
function bin_events(::FITSData, events, channel_range::NTuple{2,<:Integer}, x_edges, y_edges)::Array{Int64}
    min_chan, max_chan = extrema(channel_range)
    @mpidebug "Binning events" extrema(channel_range) extrema(x_edges) extrema(y_edges)

    events = events[minimum(x_edges).<events[:, 1].<maximum(x_edges), :]
    events = events[minimum(y_edges).<events[:, 2].<maximum(y_edges), :]
    events = events[min_chan.<events[:, 3].<max_chan, :]

    # Maximum doesn't like the Unitful array
    channels = trunc.(Int64, ustrip.(events[:, 3]))

    n_channels = max_chan - min_chan + 1 # +1 because last channel is counted
    @mpidebug "Channel data" extrema(channels)

    binned = zeros(Int64, n_channels, length(x_edges) - 1, length(y_edges) - 1)

    MPI.Init()
    if MPI.Comm_rank(comm) == 0
        is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
        prog = ProgressUnknown("Events read:", enabled=!is_logging(stderr), dt=1)
    end
    for r in eachrow(events)
        i = searchsortedlast(x_edges, r[1])
        j = searchsortedlast(y_edges, r[2])
        # shift channel number
        chan = trunc(Int64, r[3]) - min_chan + 1 # +1 because we want count to start at 1
        binned[chan, i, j] += 1
        if MPI.Comm_rank(comm) == 0
            ProgressMeter.next!(prog)
        end
    end
    if MPI.Comm_rank(comm) == 0
        ProgressMeter.finish!(prog)
    end

    return binned
end