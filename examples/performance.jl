using Profile, BenchmarkTools
using Unitful, UnitfulAstro
using BayesJ

println("PROFILE GENERATION")
nfw_params = (5.0e14, 0.13, 4.0, 1.0620, 5.4807, 0.3292, 1.156)
einasto_params = (5.0e14, 0.13, 4.0, 1.0, 1.0620, 5.4807, 0.3292, 1.156)
v2006_params = (4.705e-3, 0.247e-1, 94.6, 75.83, 0.916, 0.526, 3.607, 4.943, 1239.9, 3.61, 0.27, 57, 3.88, 1.42, 0.12, 5.00, 10.0)

println("Model_NFW")
display(@benchmark Model_NFW($nfw_params..., z=0.1))

println("\nModel_Einasto")
display(@benchmark Model_Einasto($einasto_params..., z=0.1))

println("\nModel_Vikhlinin2006")
display(@benchmark Model_Vikhlinin2006($v2006_params...))


println("\n\nPOINT CALCULATIONS")
nfw_temperature, nfw_density = Model_NFW(nfw_params..., z=0.1)
einasto_temperature, einasto_density = Model_Einasto(einasto_params..., z=0.1)
v2006_temperature, v2006_density = Model_Vikhlinin2006(v2006_params...)

println("Model_NFW")
println("Temperature")
display(@benchmark nfw_temperature($rand() * 1u"Mpc"))
println("Density")
display(@benchmark nfw_density($rand() * 1u"Mpc"))

println("\nModel_Einasto")
println("Temperature")
display(@benchmark einasto_temperature($rand() * 1u"Mpc"))
println("Density")
display(@benchmark einasto_density($rand() * 1u"Mpc"))

println("\nModel_Vikhlinin2006")
println("Temperature")
display(@benchmark v2006_temperature($rand() * 1u"Mpc"))
println("Density")
display(@benchmark v2006_density($rand() * 1u"Mpc"))