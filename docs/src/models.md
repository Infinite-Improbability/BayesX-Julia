# [Cluster Models](@id cluster_models_page)

When setting up sampling the user must supply a cluster model function to the [`sample`](@ref) function. This model takes a set of parameters (typically drawn from the priors) and returns two functions for the gas temperature and gas mass density as a function of radius.

## Included Models

BayesJ has three cluster models included. These are the [`Model_NFW`](@ref), [`Model_Einasto`](@ref) and [`Model_Vikhlinin2006`](@ref).

### NFW

The NFW model [`Model_NFW`](@ref) assumes the dark matter density profile follows the Navarro–Frenk–White profile and the electron pressure profile follows the generalised NFW model.

This model is derived in [olamaieSimpleParametricModel2012](@cite). In additional to the aforementioned profiles it assumes
- The cluster is spherical
- The ICM is in hydrostatic equilbrium
- The ICM is an ideal gas
- The local gas fraction (gas mass / total mass) is much less than unity at all radii.

#### Parameters
The following parameters can be investigated with priors.

| Parameter  | Definition                                           |
|------------|------------------------------------------------------|
| MT_Δ     | Total mass within the overdensity radius ``R_{Δ}`` |
| fg_Δ     | Gas fraction at ``R_{Δ}``                          |
| c_Δ_dm   | NFW concentration parameter                          |
| α          | GNFW parameter                                       |
| β          | GNFW parameter                                       |
| γ          | GNFW parameter                                       |
| c_Δ_GNFW | GNFW gas concentration parameter                     |

### Einasto

The Einasto model [`Model_Einasto`](@ref) assumes the dark matter density profile follows the Einasto profile and the gas pressure profile follows the generalised NFW model.

This model makes the same assumptions as the [NFW](@ref) model but subsitutes the Einasto profile in place of the NFW profile for dark matter density.

This model is derived in [olamaieBAYESXBayesianInference2015](@cite).

#### Parameters
The following parameters can be investigated with priors.

| Parameter    | Definition                                           |
|--------------|------------------------------------------------------|
| MT_Δ       | Total mass within the overdensity radius ``R_{Δ}`` |
| fg_Δ       | Gas fraction at ``R_{Δ}``                          |
| c_Δ_dm     | NFW concentration parameter                          |
| n            | Einasto index                                        |
| α            | GNFW parameter                                       |
| β            | GNFW parameter                                       |
| γ            | GNFW parameter                                       |
| c\_Δ\_GNFW | GNFW gas concentration parameter                     |

### Vikhlinin 2006

The [`Model_Vikhlinin2006`](@ref) is the model used by Vikhlinin et al. in an analysis of Chandra observations [vikhlininChandraSampleNearby2006](@cite).

This model is highly flexible, featuring independent equations for gas
number density profile and gas temperature profile.

#### Parameters
The following parameters can be investigated with priors.

| Parameter | Equation     | Definition                                |
|-----------|--------------|-------------------------------------------|
| n0        | Density      | Number density normalisation              |
| n02       | Density      | Number density normalisation of core term |
| α         | Density      | Dimensionless parameter                   |
| β         | Density      | Dimensionless parameter                   |
| β2        | Density      | Dimensionless parameter                   |
| ϵ         | Density      | Dimensionless parameter                   |
| rc        | Density      | Scale radius                              |
| rc2       | Density      | Scale radius of core term                 |
| rs        | Density      | Scale radius                              |
| T0        | Temperature  | Temperature normalisation                 |
| TminT0    | Temperature  | Ratio of Tmin to T0                       |
| rcool     | Temperature  | Scale radius                              |
| acool     | Temperature  | Dimensionless parameter                   |
| rt        | Temperature  | Scale radius                              |
| a         | Temperature  | Dimensionless parameter                   |
| b         | Temperature  | Dimensionless parameter                   |
| c         | Temperature  | Dimensionless parameter                   |









## Custom Models

The user can easily subsitute a custom model. It should be a function that will take a sequence of parameters associated with the priors as positional arguements and accept the keyword argument `z` for redshift. This value can then be ignored if it is of no consequence to the model as is the case in [`Model_Vikhlinin2006`](@ref). As UltraNest cannnot accept priors with units a version of the function may need to be provided to handle cases with/without units.

The function should return two functions, for gas temperature and gas mass density. Each function should take a single arguement: radius (on the sky) with respect to the cluster center. They should return a single scalar value for the associated property. Both the radius and output should have appropriate units using the system implemented by `Unitful.jl`. The temperature should be in units of energy  (so more strictly it is ``k_B T``).

The [`priorcheck`](@ref) function is provided to help with applying constraints to input priors during sampling.