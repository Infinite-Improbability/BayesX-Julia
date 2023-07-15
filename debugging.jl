using Plots
using Profile
using BenchmarkTools

include("src/run.jl")

data = FITSData(
    "/home/ryan/data/chandra/4361/manual3/repro/acisf04361_repro_evt2.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/bg_trimmed_300-7000.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.arf",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.rmf",
    0.492u"arcsecond"
)

response = load_response(data, [0.3u"keV", 7u"keV"])

em = prepare_model_mekal(3.89, 0.1, (0.3:0.01:7)u"keV")

model = Model_NFW_GNFW(1e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.164, [24, 24], 0.492u"arcsecond", em, 100e3u"s", response)
@profview Model_NFW_GNFW(1e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.164, [24, 24], 0.492u"arcsecond", em, 100e3u"s", response)
display(@benchmark Model_NFW_GNFW(1e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.164, [24, 24], 0.492u"arcsecond", em, 100e3u"s", response))
# @code_warntype Model_NFW_GNFW(1e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.164, [24, 24], 0.492u"arcsecond", em, 100e3u"s", response)

hm(i) = heatmap(model[i, :, :])

display(hm(1))

display(heatmap(dropdims(sum(model, dims=1), dims=1)))