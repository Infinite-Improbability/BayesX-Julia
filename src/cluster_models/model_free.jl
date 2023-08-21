export Model_Free

function Model_Free(params...; z)
    @argcheck length(params) % 2 == 0

    overdensities = [200, 500, 2500, 1000, 100]
    @argcheck length(params) <= length(overdensities)

    radii = [r * 1u"kpc" for r in params[1:2:end]]
    pressures = [p * 1u"Pa" for p in params[2:2:end]]

    display(radii)
    display(pressures)

    pressure = extrapolate(Interpolations.scale(interpolate(pressures, BSpline(Quadratic(Reflect(OnCell())))), radii), Line())

    ρ_crit_z = ρ_crit(z)
    mass(overdensity, r) = mass = 4π / 3 * overdensity * ρ_crit_z * r^3

    function density_at_overdensity(overdensity, r)
        dPdr = Interpolations.gradient(pressure, r)
        m = mass(overdensity, r)
        -dPdr * r^2 * G * m # hydrostatic equilibrium
    end

    densities = [density_at_overdensity(o, r) for (o, r) in zip(overdensities, radii)]
    density = extrapolate(Interpolations.scale(interpolate(densities, BSpline(Quadratic(Reflect(OnCell())))), radii), Line())

    temperature(r) = μ * pressure(r) / k_B / density(r) # ideal gas

    return temperature, pressure
end