using SpectralFitting

@xspecmodel :C_mekal struct XS_Mekal{T,F} <: SpectralFitting.AbstractSpectralModel{T,SpectralFitting.Additive}
    "Normalisation"
    K::T
    "Plasma Temperature"
    t::T
    "Hydrogen Density"
    ρ::T
    "Metal Abundances"
    a::T
    "Redshift"
    z::T
    "Switch"
    s::Int
end
function XS_Mekal(;
    K=FitParam(1.0),
    t=FitParam(8.0),
    ρ=FitParam(10.0),
    a=FitParam(1.0),
    z=FitParam(0.1),
    s=0
)
    XS_Mekal{typeof(K),SpectralFitting.FreeParameters{(:K, :t, :ρ)}}(
        K, t, ρ, a, z, s
    )
end
SpectralFitting.register_model_data(XS_Mekal, "mekal1.dat", "mekal2.dat", "mekal3.dat", "mekal4.dat", "mekal5.dat", "mekal6.dat")


# mekal = XS_Mekal(t=FitParam(8.0), ρ=FitParam(12.0), z=FitParam(0.1))
# invokemodel(collect(0.1:0.1:2), mekal)
# invokemodel(collect(0.1:0.1:2), mekal, [1.0, 6.0, 4.0])

function prepare_model_mekal(
    nHcol,
    redshift,
    energy_bins;
    temperatures=1e-30:0.01:15,
    densities=0:10.0:100,
    normalisation=1.0
)
    @info "Preparing model"

    energy_bins = collect(energy_bins)

    # TODO: Figure out normalisation

    # Generate transmission fractions
    absorption_model = PhotoelectricAbsorption(FitParam(nHcol))
    absorption = invokemodel(energy_bins, absorption_model)

    @assert all(isfinite, absorption)

    # Generate source flux
    # TODO: unit?!
    emission_model(kbT, ρ) = XS_Mekal(K=FitParam(normalisation), t=FitParam(kbT), ρ=FitParam(ρ), z=FitParam(redshift))
    points = [emission_model(t, d) for t in temperatures, d in densities]
    emission = invokemodel.(
        Ref(energy_bins),
        points
    )

    @assert all(all.(isfinite, emission))

    # Apply absorption
    flux = [absorption .* emission[index] for index in eachindex(IndexCartesian(), emission)]

    @assert all(all.(isfinite, flux))

    return linear_interpolation(
        (temperatures, densities),
        flux
    )
end
# const surrogate = prepare_model(2.2, 0.1, 0.3:0.1:3.0)
