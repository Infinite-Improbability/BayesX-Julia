using ArgCheck
using Unitful, UnitfulAstro, DimensionfulAngles
using Integrals

include("params.jl")
include("emission.jl")

"""
    p_crit(z)

Calculate the critical density at some redshift `z`.
"""
ρ_crit(z) = 3 * H(cosmo, z)^2 / (8π * G)

"""
    PriorError(likelihood)

This wraps a fallback likelihood value so it can be passed up the chain on an invalid
prior combination.
"""
struct PriorError <: Exception
    likelihood::Float64
end

"""
    priorcheck(condition, likelihood)

Check that a condition holds and fallback to the given
likelihood if not.

The likelihood is encapsulated in a [PriorError](@ref) so a try-catch
block in the function calling the model can extract it and
return it as the likelihood without going throw observation generation
and likelihood calculation.

The faux likelihood should be chosen to be sufficiently small (e.g. -1e100) that it
will not be competing with the likelihood from good parameters. It should also be
designed with a slope to point the sampler in the right direction. For instance if 
the constraint is that ``a < b`` then ``-1e100 * (1 + (a - b))`` would  be good.
This constraint ensures the likelihood grows larger when `a` and `b` grow closer 
together. The additional one ensures that ``a=b`` won't set the likelihood to zero.


See the [ultranest documentation](https://johannesbuchner.github.io/UltraNest/priors.html#Complicated-constraints-and-rejection-in-the-likelihood) 
for additional details.
"""
function priorcheck(condition, likelihood)
    if !condition
        throw(PriorError(likelihood))
    end
end

include("cluster_models/model_nfw.jl")
include("cluster_models/model_einasto.jl")
include("cluster_models/model_vikhlinin2006.jl")

