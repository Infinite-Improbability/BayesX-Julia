# Cluster Models

When setting up sampling the user must supply a cluster model function to the [sample](@ref) function.
This model takes a set of parameters and returns two functions for the gas temperature and gas mass density as a function of
radius.

## Included Models

BayesJ has three cluster models included. These are the [Model_NFW](@ref), [Model_Einasto](@ref) and [Model_Vikhlinin2006](@ref).

### NFW

The NFW model [Model_NFW](@ref) assumes the dark matter density profile follows the Navarro-Frenk-White profile and the gas pressure profile follows the
generalised NFW model.



### Einasto

The Einasto model [Model_Einasto](@ref) assumes the dark matter density profile follows the Einasto profile and the gas pressure profile follows the
generalised NFW model.

### Vikhlinin 2006

The [Model_Vikhlinin2006](@ref) is based on a model used by Vikhlinin et al. in an analysis of Chandra observations.

TODO: Details, citations for all models. Advantages and disadvantages

## Custom Models

The user can easily subsitute a custom model. It should be a function that will take a sequence of parameters associated with the priors as positional arguements and accept the keyword argument `z` for redshift. This value can then be ignored if it is of no consequence to the model as is the case in `Vikhlinin2006``. As UltraNest cannnot accept priors with units a version of the function may need to be provided to handle cases with/without units.

The function should return two functions, for gas temperature and gas mass density. Each function should take a single arguement: radius (on the sky) with respect to the cluster center. They should return a single scalar value for the associated property. Both the radius and output should have appropriate units using the system implemented by `Unitful.jl`. The temperature should be in units of energy  (so more strictly it is ``k_B T``).

The [priorcheck](@ref) function is provided to help with applying
constraints to input priors during sampling.