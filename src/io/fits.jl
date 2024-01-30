using FITSIO

export FITSData

# Import our patch
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
    @mpiwarn "You should use dimensionful angles"
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
Assumes the RMF and ARF are formatted according to the HEASARC OGIP Calibration Memo CAL/GEN/92-002
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

    # Get the HDU with the matrix
    rmf_hdus::Vector{TableHDU} = [
        h for h in f_rmf if safe_read_key(h, "extname", "Exception when looking for matrix HDU. Probably a HDU without extname.")[1] in ["MATRIX", "SPECRESP MATRIX"]
    ]
    if length(rmf_hdus) > 1
        @mpiwarn "$(length(rmf_hdus)) HDUs with matrix extension found. Using the first."
    end
    rmf = rmf_hdus[1]
    @mpidebug "Selected HDU with MATRIX extension and HDUNAME '$(safe_read_key(rmf, "HDUNAME", "HDU has no name")[1])' from $(data.rmf) for RMF"

    # Load the header
    header = read_header(rmf)

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
    # Both are given as vectors because there can be multiple ranges
    # This is to optimise the storage of sparse matrices
    # And also to annoy me
    # (See https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html)
    # Furthermore they may use fixed or variable length arrays

    # We need to find the index of the F_CHAN column
    # It should always be 4 but we won't assume
    f_index = 0
    for key in keys(header)
        if header[key] == "F_CHAN"
            stripped = strip(key, ('T', 'Y', 'P', 'E'))
            f_index = parse(Int, stripped)
            break
        end
    end
    if f_index == 0
        @mpierror "Could not find F_CHAN column in RMF header"
    elseif f_index != 4
        @mpiwarn "F_CHAN column is not in the expected position (4). This suggests a nonstandard RMF."
    end

    n_index = 0
    for key in keys(header)
        if header[key] == "N_CHAN"
            stripped = strip(key, ('T', 'Y', 'P', 'E'))
            n_index = parse(Int, stripped)
            break
        end
    end
    if n_index == 0
        @mpierror "Could not find N_CHAN column in RMF header"
    elseif n_index != 5
        @mpiwarn "N_CHAN column is not in the expected position (5). This suggests a nonstandard RMF."
    end

    m_index = 0
    for key in keys(header)
        if header[key] == "MATRIX"
            stripped = strip(key, ('T', 'Y', 'P', 'E'))
            m_index = parse(Int, stripped)
            break
        end
    end
    if m_index == 0
        @mpierror "Could not find MATRIX column in RMF header"
    elseif m_index != 6
        @mpiwarn "MATRIX column is not in the expected position (6). This suggests a nonstandard RMF."
    end

    # Load in the channel ranges
    first_channels = read(rmf, "F_CHAN")
    num_channels = read(rmf, "N_CHAN")

    # But fixed and variable length arrays behave differently
    # We add special fixed length handling here
    # Format code examples follow
    # Fixed length: 3J
    # Variable length: PJ(3)
    f_format = header["TFORM$f_index"] # First channel data format
    if isdigit(f_format[1])
        digits = []
        for c in f_format
            if isdigit(c)
                push!(digits, c)
            else
                break
            end
        end
        len = parse(Int, join(digits))

        if 'P' in f_format && len > 1
            @mpierror "F_CHAN reports a fixed length array of length $len but also has a variable length indicator, P, in the format. This is nonstandard and not supported."
        end

        first_channels = [collect(Iterators.partition(first_channels, len))]
    end

    n_format = header["TFORM$n_index"]
    if isdigit(n_format[1])
        digits = []
        for c in n_format
            if isdigit(c)
                push!(digits, c)
            else
                break
            end
        end
        len = parse(Int, join(digits))

        if 'P' in n_format && len > 1
            @mpierror "N_CHAN reports a fixed length array of length $len but also has a variable length indicator, P, in the format. This is nonstandard and not supported."
        end

        num_channels = [collect(Iterators.partition(num_channels, len))]
    end

    @assert all(length.(first_channels) .== length.(num_channels))
    # Calculate the last channel as first_channel + number_of_channels - 1
    # The -1 accounts for the first channel channel contributing to the total number of channels
    last_channels = [f .+ n .- 1 for (f, n) in zip(first_channels, num_channels)]

    # Preallocate output matrix
    # Response matrix has shape (PI, E) where PI is channel and E is energy bin
    # NAXIS gives the number of rows (energy bins) in the table
    # We'll verify this matches our other data
    n_bins = header["NAXIS2"]
    @assert n_bins == length(first_channels) "Claimed number of rows in table ($n_bins) does not match length of data $(length(first_channels))"
    # Get the number of channels from the header
    n_channels = header["DETCHANS"]
    response_matrix = zeros(Float64, n_channels, n_bins)

    # Then we load in MATRIX
    # For each bin this a vector giving the probability of a photon going to each channel
    # So we want the vector to be the same length as the number of channels in the bin
    # How does this interact with multiple range? Presumbly it should be the sum of the channel subsets.
    energy_to_channel_mapping = read(rmf, "MATRIX")
    m_format = header["TFORM$m_index"]
    if isdigit(m_format[1])
        digits = []
        for c in m_format
            if isdigit(c)
                push!(digits, c)
            else
                break
            end
        end
        len = parse(Int, join(digits))

        if 'P' in m_format && len > 1
            @mpierror "MATRIX reports a fixed length array of length $len but also has a variable length indicator, P, in the format. This is nonstandard and not supported."
        end

        energy_to_channel_mapping = [collect(Iterators.partition(energy_to_channel_mapping, len))]
    end

    @assert length.(energy_to_channel_mapping) == sum.(num_channels) "Discrepency between stated number of channels in bin (N_CHAN) and length of energy to channel mapping vector (MATRIX)"

    e_to_c = Vector{Vector{Vector{Float64}}}(undef, n_bins)
    for i in eachindex(energy_to_channel_mapping)
        ns = num_channels[i]
        ec = energy_to_channel_mapping[i]
        pointer = 1
        e_to_c[i] = Vector{Vector{Float64}}(undef, length(ns))
        for j in 1:length(ns) # we don't want eachindex(ns) because we want order guaranteed
            n = ns[j]
            last = pointer + n - 1
            e_to_c[i][j] = ec[pointer:last]
            pointer = last + 1
        end
        @assert pointer == length(ec) + 1
    end

    # tidy up
    close(f_rmf)

    # Load values into matrix
    # We iterate over energy bins, which are our columns
    @views for bin in axes(response_matrix, 2)
        # Each bin may have several ranges
        for (f, l, ec) in zip(first_channels[bin], last_channels[bin], e_to_c[bin])
            # for each bin `i` we select the appropriate channel range from the matrix with first_channel_of_bin[i]:last_channel_of_bin[i]
            # and assign the energy mapping to it
            response_matrix[f:l, bin] = ec
        end
    end

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

    # Get valid bin range
    # Get first bin where lower edge >= user specified `min_energy``
    min_bin = searchsortedfirst(rmf_energy_bin_minimums, min_energy)
    # Get last bin where upper edge <= user specificed `max_energy``
    max_bin = searchsortedlast(rmf_energy_bin_maximums, max_energy)

    # Get minimum and maximum channels in valid bin range
    min_channel = minimum(minimum(first_channels[min_bin:max_bin]))
    max_channel = maximum(maximum(last_channels[min_bin:max_bin]))

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