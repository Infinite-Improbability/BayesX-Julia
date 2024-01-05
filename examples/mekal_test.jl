using Plots
using Profile
using BenchmarkTools
using Unitful, DimensionfulAngles, UnitfulAstro
using PoissonRandom

# include("../src/BayesJ.jl")

using BayesJ

# ENV["JULIA_DEBUG"] = "BayesJ"

energy_range = range(0.3u"keV", 7.0u"keV", 33)
redshift = 0.164

em = BayesJ.prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_range,
    redshift,
    use_interpolation=false,
)

T = 0.0u"keV"
ρ = 0.0u"g/cm^3"
ΔT = 1e-3u"keV"
Δρ = 1e-3u"g/cm^3"
loops = 0

while loops < 10
    global loops
    global ΔT
    global Δρ
    global T
    global ρ

    loops = 0
    T = 1.0e-10u"keV"
    # ρ = 4.4e-5u"g/cm^3"
    ρ = 0.0u"g/cm^3"

    # ΔT /= 10
    Δρ /= 10

    counts = em(T, BayesJ.hydrogen_number_density(ρ))

    while iszero(counts)
        global loops
        global T
        global ρ
        # T += ΔT
        ρ += Δρ
        counts .= em(T, BayesJ.hydrogen_number_density(ρ))
        loops += 1
        # if loops % 100 == 0
        #     @info "Loop progressing" loops T ρ ΔT Δρ
        # end
    end

    @info "Inner loop ended" loops T ρ ΔT Δρ
end

@info "Found first nonzero value at" loops T ρ uconvert(u"cm^-3", BayesJ.hydrogen_number_density(ρ))