using Plots
using Profile
using BenchmarkTools
# using FortranFiles
using Unitful, DimensionfulAngles

include("src/run.jl")

data = FITSData(
    "/home/ryan/data/chandra/4361/manual3/repro/acisf04361_repro_evt2.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/bg_trimmed_300-7000.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.arf",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.rmf",
    0.492u"arcsecondᵃ"
)

energy_range = (0.3:0.01:0.46)u"keV"
response = load_response(data, energy_range)
em = prepare_model_mekal2(3.89, energy_range)
model = Model_NFW_GNFW(5e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.5, [32, 32], 0.492u"arcsecondᵃ", em, 300e3u"s", response)

# f = FortranFile("../BayesX/log.dat")
# fm = reshape(read(f, (Float64, 64 * 64 * 32)), (32, 64, 64))