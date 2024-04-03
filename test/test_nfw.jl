using BayesJ
using Test
using Unitful, UnitfulAstro
using BayesJ.ModelNFW
import Cosmology.a2E

function test()
    mass = 5e14u"Msun"
    fg = 0.13
    c_dm = 3.0
    gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
    z = 1
    r = 1u"Mpc"

    @testset "Helper Functions" begin
        ρ_crit_z = ModelNFW.ρ_crit(z)
        r_Δ = ModelNFW.r_delta(mass, 500, ρ_crit_z)
        r_s = uconvert(u"Mpc", r_Δ / c_dm)
        ρ_s = ModelNFW.rho_s(ρ_crit_z, 500, c_dm)
        r_p = uconvert(u"Mpc", r_Δ/gnfw[4])
        pei_integrand = ModelNFW.gnfw_gas_mass_integrand(r, r_s, r_p, gnfw[1], gnfw[2], gnfw[3])
        @test ρ_crit_z ≈ 25963120615613.40u"kg*km^2/m^3/Mpc^2" rtol=1e-12
        @test r_Δ ≈ 8.39917483737E-01u"Mpc" rtol=1e-12
        @test ρ_s ≈ 1.83616341601E+17u"kg*km^2/m^3/Mpc^2" rtol=1e-12
        @test r_p ≈ 7.13608737245E-01u"Mpc" rtol=1e-12
        @test pei_integrand ≈  5.01482647641E-01 broken=true
    end
    
end

test()