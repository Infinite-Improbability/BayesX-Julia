using BayesJ
using BenchmarkTools
using Unitful

energy_bins = range(0u"keV", 13u"keV", length=1000)

emission_model = BayesJ.prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_bins,
    0.1,
    use_interpolation=false
)

T = 1.0u"keV"
nH = 1.0e-3u"cm^-3"

display(@benchmark emission_model($T, $nH))