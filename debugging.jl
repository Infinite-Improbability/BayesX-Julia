using Plots

include("src/run.jl")

em = prepare_model_mekal(3.89, 0.1, (0.3:0.01:7)u"keV")

model = Model_NFW_GNFW(1e14u"Msun", 0.13, 1.062, 5.4807, 0.3292, 1.156, 0.164, [24, 24], 0.492u"arcsecond", em, 100e3u"s", [1 0; 0 1])

hm(i) = heatmap([p[i] for p in model])

hm(1)

heatmap([sum(p) for p in model])