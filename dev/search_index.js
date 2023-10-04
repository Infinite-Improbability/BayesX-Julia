var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"","category":"page"},{"location":"private/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"## Index\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@index Pages = [\"private.md\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n## Private Interface\n\n### Cluster Models\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"gas_models.jl\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n### Constants\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"params.jl\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n### Emissions Modelling\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"emission.jl\", \"mekal.jl\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n### IO\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"io.jl\", \"fitsio_fix.jl\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n### Likelihood\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"likelihood.jl\"]","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"\n### Utilities\n","category":"page"},{"location":"private/","page":"Internals","title":"Internals","text":"@autodocs Modules = [BayesJ] Public  = false Private = true Pages   = [\"mpi.jl\", \"run.jl\"] ```","category":"page"},{"location":"models/#Cluster-Models","page":"Cluster Models","title":"Cluster Models","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"When setting up sampling the user must supply a cluster model function to the sample function. This model takes a set of parameters (typically drawn from the priors) and returns two functions for the gas temperature and gas mass density as a function of radius.","category":"page"},{"location":"models/#Included-Models","page":"Cluster Models","title":"Included Models","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"BayesJ has three cluster models included. These are the Model_NFW, Model_Einasto and Model_Vikhlinin2006.","category":"page"},{"location":"models/#NFW","page":"Cluster Models","title":"NFW","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The NFW model Model_NFW assumes the dark matter density profile follows the Navarro–Frenk–White profile and the electron pressure profile follows the generalised NFW model.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"This model is derived in (Olamaie et al., 2012). In additional to the aforementioned profiles it assumes","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The cluster is spherical\nThe ICM is in hydrostatic equilbrium\nThe ICM is an ideal gas\nThe local gas fraction (gas mass / total mass) is much less than unity at all radii.","category":"page"},{"location":"models/#Parameters","page":"Cluster Models","title":"Parameters","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The following parameters can be investigated with priors.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"Parameter Definition\nMT_Δ Total mass within the overdensity radius R_Δ\nfg_Δ Gas fraction at R_Δ\ncΔdm NFW concentration parameter\nα GNFW parameter\nβ GNFW parameter\nγ GNFW parameter\ncΔGNFW GNFW gas concentration parameter","category":"page"},{"location":"models/#Einasto","page":"Cluster Models","title":"Einasto","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The Einasto model Model_Einasto assumes the dark matter density profile follows the Einasto profile and the gas pressure profile follows the generalised NFW model.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"This model makes the same assumptions as the NFW model but subsitutes the Einasto profile in place of the NFW profile for dark matter density.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"This model is derived in (Olamaie et al., 2015).","category":"page"},{"location":"models/#Parameters-2","page":"Cluster Models","title":"Parameters","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The following parameters can be investigated with priors.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"Parameter Definition\nMT_Δ Total mass within the overdensity radius R_Δ\nfg_Δ Gas fraction at R_Δ\ncΔdm NFW concentration parameter\nn Einasto index\nα GNFW parameter\nβ GNFW parameter\nγ GNFW parameter\nc_Δ_GNFW GNFW gas concentration parameter","category":"page"},{"location":"models/#Vikhlinin-2006","page":"Cluster Models","title":"Vikhlinin 2006","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The Model_Vikhlinin2006 is the model used by Vikhlinin et al. in an analysis of Chandra observations (Vikhlinin et al., 2006).","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"This model is highly flexible, featuring independent equations for gas number density profile and gas temperature profile.","category":"page"},{"location":"models/#Parameters-3","page":"Cluster Models","title":"Parameters","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The following parameters can be investigated with priors.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"Parameter Equation Definition\nn0 Density Number density normalisation\nn02 Density Number density normalisation of core term\nα Density Dimensionless parameter\nβ Density Dimensionless parameter\nβ2 Density Dimensionless parameter\nϵ Density Dimensionless parameter\nrc Density Scale radius\nrc2 Density Scale radius of core term\nrs Density Scale radius\nT0 Temperature Temperature normalisation\nTminT0 Temperature Ratio of Tmin to T0\nrcool Temperature Scale radius\nacool Temperature Dimensionless parameter\nrt Temperature Scale radius\na Temperature Dimensionless parameter\nb Temperature Dimensionless parameter\nc Temperature Dimensionless parameter","category":"page"},{"location":"models/#Custom-Models","page":"Cluster Models","title":"Custom Models","text":"","category":"section"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The user can easily subsitute a custom model. It should be a function that will take a sequence of parameters associated with the priors as positional arguements and accept the keyword argument z for redshift. This value can then be ignored if it is of no consequence to the model as is the case in Model_Vikhlinin2006. As UltraNest cannnot accept priors with units a version of the function may need to be provided to handle cases with/without units.","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The function should return two functions, for gas temperature and gas mass density. Each function should take a single arguement: radius (on the sky) with respect to the cluster center. They should return a single scalar value for the associated property. Both the radius and output should have appropriate units using the system implemented by Unitful.jl. The temperature should be in units of energy  (so more strictly it is k_B T).","category":"page"},{"location":"models/","page":"Cluster Models","title":"Cluster Models","text":"The priorcheck function is provided to help with applying constraints to input priors during sampling.","category":"page"},{"location":"public/#Public-Documentation","page":"Public","title":"Public Documentation","text":"","category":"section"},{"location":"public/#Index","page":"Public","title":"Index","text":"","category":"section"},{"location":"public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Public-Interface","page":"Public","title":"Public Interface","text":"","category":"section"},{"location":"public/","page":"Public","title":"Public","text":"Modules = [BayesJ]\nPublic  = true\nPrivate = false","category":"page"},{"location":"public/#BayesJ.DeltaPrior","page":"Public","title":"BayesJ.DeltaPrior","text":"DeltaPrior(name::AbstractString, value::Number)\n\nA delta prior that always returns a constant value.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.FITSData","page":"Public","title":"BayesJ.FITSData","text":"FITSData(observation::AbstractString, background::AbstractString, arf::AbstractString, rmf::AbstractString, pixel_edge_angle::DimensionfulAngles.Angle)\n\nStore a collection of observation data from fits files.\n\nTakes paths to the files as inputs.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.GenericPrior","page":"Public","title":"BayesJ.GenericPrior","text":"GenericPrior(name::AbstractString, dist::Distributions.UnivariateDistribution)\n\nA generic prior that draws from the specified distribution.\n\nSee Distributions.jl for a list of avaliable distributions.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.LogUniformPrior","page":"Public","title":"BayesJ.LogUniformPrior","text":"LogUniformPrior(name::AbstractString, min::Number, max::Number)\n\nA log uniform prior that draws from a distribution between min and max whose base 10 logarithm is uniformly distributed.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.NormalPrior","page":"Public","title":"BayesJ.NormalPrior","text":"NormalPrior(name::AbstractString, mean::Number, σ::Number)\n\nA normal prior that draws from a Gaussian/normal distribution with mean and standard deviation σ.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.PlaintextData","page":"Public","title":"BayesJ.PlaintextData","text":"PlaintextData(observation::AbstractString, background::AbstractString, arf::AbstractString, rmf::AbstractString, shape::NTuple{3,Int}, rmf_shape::NTuple{2,Int}, obs_exposure_time::Unitful.Time, bg_exposure_time::Unitful.Time, pixel_edge_angle::DimensionfulAngles.Angle, rmf_unit::Unitful.Area)\n\nStore preprocessed data from plaintext files with whitespace delimiters.\n\nEvents and response files provided as paths. This data is assumed to have the appropriate binning and energy ranges already applied. The shape of the event arrays should be (channels, x, y) and the rmf_shape should be (channels, bins).\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.PriorError","page":"Public","title":"BayesJ.PriorError","text":"PriorError(likelihood)\n\nThis wraps a fallback likelihood value so it can be passed up the chain on an invalid prior combination.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.UniformPrior","page":"Public","title":"BayesJ.UniformPrior","text":"UniformPrior(name::AbstractString, min::Number, max::Number)\n\nA uniform prior that draws from a uniform distribution between min and max.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.Model_Einasto-Tuple{Real, Vararg{Any, 7}}","page":"Public","title":"BayesJ.Model_Einasto","text":"Model_Einasto(MT_Δ::Real, fg_Δ, n, α, β, γ, c_Δ_GNFW; z, Δ=500)\n\nUnitless wrapper for Model_Einasto. Mass is in solar masses.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_Einasto-Tuple{Union{Unitful.Quantity{T, 𝐌, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐌, U}} where {L, S}} where {T, U}, Vararg{Any, 7}}","page":"Public","title":"BayesJ.Model_Einasto","text":"Model_Einasto(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, n, α, β, γ, c_Δ_GNFW; z, Δ=500)\n\nGenerate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.\n\nSee Appendix A of (Olamaie et al., 2015) for a derivation of the model based on these profiles.\n\nReturns functions for gas temperature and gas mass density as a function of radius.\n\nIf n>2 then the gas temperature starts increasing to physically improbable levels (orders exceeding 10¹¹ keV). We thus constrain the value to below that. \n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_NFW-Tuple{Real, Vararg{Any, 6}}","page":"Public","title":"BayesJ.Model_NFW","text":"Model_NFW(MT_Δ::Real, fg_Δ, α, β, γ, c_Δ_GNFW; z, Δ=500)\n\nUnitless wrapper for Model_NFW. Mass is in solar masses.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_NFW-Tuple{Union{Unitful.Quantity{T, 𝐌, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐌, U}} where {L, S}} where {T, U}, Vararg{Any, 6}}","page":"Public","title":"BayesJ.Model_NFW","text":"Model_NFW(MT_Δ::Unitful.Mass, fg_Δ, c_Δ_dm, α, β, γ, c_Δ_GNFW; z, Δ=500)\n\nGenerate a cluster profile based on the NFW mass density and GNFW gas density profiles.\n\nUses the model from (Olamaie et al., 2012), which is based on the NFW dark matter density profile and the GNFW gas pressure profile.\n\nKeyword arguments are redshift and overdensity.\n\nReturns functions for gas temperature and gas mass density as a function of radius.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_Vikhlinin2006-Tuple{Real, Real, Real, Real, Any, Any, Any, Any, Real, Real, Any, Real, Any, Real, Any, Any, Any}","page":"Public","title":"BayesJ.Model_Vikhlinin2006","text":"Model_Vikhlinin2006(n0::Real, n02::Real, rc::Real, rc2::Real, α, β, β2, ϵ, rs::Real, T0::Real, TminT0, rcool::Real, acool, rt::Real, a, b, c; γ=3, kwargs...)\n\nUnitless wrapper for Model_Vikhlinin2006. Densities are in cm^-3 and radii are in kpc. Temperature is in keV.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_Vikhlinin2006-Tuple{Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Any, Any, Any, Any, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}} where {L, S}} where {T, U}, Any, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Any, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Any, Any, Any}","page":"Public","title":"BayesJ.Model_Vikhlinin2006","text":"Model_Vikhlinin2006(n0::NumberDensity, n02::NumberDensity, rc::Unitful.Length, rc2::Unitful.Length, α, β, β2, ϵ, rs::Unitful.Length, T0::Unitful.Energy, TminT0, rcool::Unitful.Length, acool, rt::Unitful.Length, a, b, c; γ=3, kwargs...)\n\nGenerate a cluster profile using the highly free models from (Vikhlinin et al., 2006).\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.priorcheck-Tuple{Any, Any}","page":"Public","title":"BayesJ.priorcheck","text":"priorcheck(condition, likelihood)\n\nCheck that a condition holds and fallback to the given likelihood if not.\n\nSimply call this function in your model. If the condition fails a PriorError is raised.\n\nThe likelihood is encapsulated in the PriorError so that a try-catch block in the function calling the model can extract it and return it as the likelihood generated by that prior set without going through observation generation and likelihood calculation.\n\nThe faux likelihood should be chosen to be sufficiently small (e.g. -1e100) that it will not be competing with the likelihood from good parameters. It should also be designed with a slope to point the sampler in the right direction. For instance if  the constraint is that a  b then -1e100 * (1 + (a - b)) would  be good. This constraint ensures the likelihood grows larger when a and b grow closer  together. The additional one ensures that a=b won't set the likelihood to zero.\n\nSee the ultranest documentation  for additional details.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.sample-Tuple{BayesJ.Dataset, AbstractRange{<:Union{Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}} where {L, S}} where {T, U}}, Function, AbstractVector{<:BayesJ.Prior}, Union{Unitful.Quantity{T, 𝐋^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-2, U}} where {L, S}} where {T, U}, Real, Tuple{var\"#s2\", var\"#s2\"} where var\"#s2\"<:Real, Tuple{var\"#s1\", var\"#s1\"} where var\"#s1\"<:Real}","page":"Public","title":"BayesJ.sample","text":"sample(data::Dataset, energy_range::AbstractRange{Unitful.Energy}, cluster_model::Function, priors::AbstractVector{Prior}, nhCol::SurfaceDensity, redshift, x, y)\n\nRun Bayesian inference on a given set of data considering only the selected energy range.\n\nAn gas emission model (density, temperature) → emissivity can be provided.\nThe first two priors should always be x0 and y0, giving cluster centre position.\nx and y are tuples of (min, max).\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.sample-Union{Tuple{T}, Tuple{T, T, Matrix, Function, Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U}, Real}} where T<:AbstractArray","page":"Public","title":"BayesJ.sample","text":"sample(observed, observed_background, response_function, transform, obs_exposure_time, bg_exposure_time, redshift; prior_names, cluster_model, emission_model, param_wrapper, pixel_edge_angle)\n\nConfigure some necessary variables and launch ultranest.\n\nThe observed array includes the background. The first dimension is energy, the other two are spatial.\nThe response function includes both the RMF and ARF, as described in apply_response_function.\nThe emission model should be a function compatible with the requirements of the surface_brightness function, which it will be passed to.\nThe pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.\n\nThis area is assumed to be square with the edge angle giving the side length.\n\nThe first two priors should always be \"x0\" and \"y0\", giving cluster centre position and the order of prior_names must match the transform.\nThe cluster model should be a function that takes the parameters (and redshift as a kwarg) and returns (gas_temperature,gas_density) as functions of\n\nradius which return their respective quantities with units.\n\nparam_wrapper takes the output of the transform function and adds any additional arguements necessary for the model.\ncentre_radius, mask and integration_limit are passed through to make_observation\n\n\n\n\n\n","category":"method"},{"location":"#BayesJ.jl-Documentation","page":"Home","title":"BayesJ.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Bayesian inference on galaxy clusters in X-ray.","category":"page"}]
}
