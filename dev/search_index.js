var documenterSearchIndex = {"docs":
[{"location":"private/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"private/#Contents","page":"Internals","title":"Contents","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Pages = [\"private.md\"]\nDepth = 3","category":"page"},{"location":"private/#Index","page":"Internals","title":"Index","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Pages = [\"private.md\"]","category":"page"},{"location":"private/#Private-Interface","page":"Internals","title":"Private Interface","text":"","category":"section"},{"location":"private/#Cluster-Models","page":"Internals","title":"Cluster Models","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"gas_models.jl\"]","category":"page"},{"location":"private/#Constants","page":"Internals","title":"Constants","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"params.jl\"]","category":"page"},{"location":"private/#BayesJ.cosmo","page":"Internals","title":"BayesJ.cosmo","text":"Cosmology assumed by BayesX\n\nDefaults to:\n\nh=0.69,\nNeff=3.04,\nOmegaK=0,\nOmegaM=0.29,\nOmegaR=nothing,\nTcmb=2.7255,\nw0=-1,\nwa=0\n\n\n\n\n\n","category":"constant"},{"location":"private/#BayesJ.μ","page":"Internals","title":"BayesJ.μ","text":"Mean mass per gas particle = 0.6mₚ\n\n\n\n\n\n","category":"constant"},{"location":"private/#BayesJ.μ_e","page":"Internals","title":"BayesJ.μ_e","text":"Mean gas mass per electron = 1.14mₚ\n\n\n\n\n\n","category":"constant"},{"location":"private/#Emissions-Modelling","page":"Internals","title":"Emissions Modelling","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"emission.jl\", \"mekal.jl\"]","category":"page"},{"location":"private/#BayesJ.apply_response_function-Union{Tuple{T}, Tuple{Vector{T}, Matrix{T}, T}} where T<:AbstractFloat","page":"Internals","title":"BayesJ.apply_response_function","text":"apply_response_function(counts_per_bin::Vector, response::Matrix, exposure_time::Unitful.Time)::Vector{Float64}\n\nApplies the response function RSP(PI E) = RMF(PI E)  ARF(E) to all energy bins and returns adjusted counts per bin.\n\nbeginaligned\nC(PI) = T int RMF(PI E)  ARF(E)  S(E)  dE \n T sum_j R_ij Aj Sj\nendaligned\n\nwith C(PI) is the observed counts in a detector channel PI, T is the observation time, ARF(E) is the effective area of the instrument and RMF(E, PI) [we use RMF(PI, E)] is the unitless response matrix. (Handbook of X-ray Astronomy Chapter 2011, Chapter 5 pg 87, by K. Arnaud, R. Smith and A. Siemiginowska)\n\nThis function takes the combined RMF and ARF as the response function. This is to recalculating it on every call. Some people format the RMF as RMF(E, PI), others as RMF(PI, E). This latter convention is used by CIAO and by us.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.make_observation-Union{Tuple{T}, Tuple{A}, Tuple{Function, Function, Any, Any, A, Any, T, Any, Tuple{Real, Real}, Any}} where {A<:(Union{Unitful.Quantity{T, 𝐀, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐀, U}} where {L, S}} where {T, U}), T<:(Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U})}","page":"Internals","title":"BayesJ.make_observation","text":"Centre position is given in arcseconds.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.make_observation-Union{Tuple{T}, Tuple{A}, Tuple{Function, Function, Any, Any, A, Any, T, Any, Tuple{var\"#s302\", var\"#s302\"} where var\"#s302\"<:(Union{Unitful.Quantity{T, 𝐀, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐀, U}} where {L, S}} where {T, U}), Any}} where {A<:(Union{Unitful.Quantity{T, 𝐀, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐀, U}} where {L, S}} where {T, U}), T<:(Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U})}","page":"Internals","title":"BayesJ.make_observation","text":"make_observation(temperature, density, z, shape, pixel_edge_angle, emission_model, exposure_time, response_function, centre, centre_radius, mask=nothing)\n\nGenerate an image of the cluster given functions for the radial profile of gas temperature and electron density and assorted observational parameters.\n\nshape describes the size of the source number as a number of spatial bins in each dimension. The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds. This area is assumed to be square with the edge angle giving the side length. The emission model should be a function compatible with the requirements of the surface_brightness function, which it will be passed to. The response function includes both the RMF and ARF, as described in apply_response_function.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.surface_brightness-Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Function, Function, Float64, Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Any, Union{Unitful.Quantity{T, 𝐀, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐀, U}} where {L, S}} where {T, U}}","page":"Internals","title":"BayesJ.surface_brightness","text":"surface_brightness(projected_radius, temperature, density, z, limit, model, pixel_edge_angle)\n\nCalculate the observed surface_brightness at some projected radius on the sky.\n\nThe temperature and density functions should take the distance from origin (cluster centre) and return the appropriate value, with units.\n\nThe model should take these functions as input and return a vector of volume emissivities, where the vector element correspond with the energy bins used to generate the mode. It should include the effects of redshift on the bins and time dilation on the count rate.\n\nThe limit controls limits on LOS radius during integration. Theoretically it should be ±infinity, but we may approximate it as a finite value.\n\nThe pixel edge angle is the angular size of a single pixel.\n\nThe result returned is the expected count rate per unit observing area, as a vector corresponding to the energy bins as used in the model.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.XS_Mekal","page":"Internals","title":"BayesJ.XS_Mekal","text":"Mekal model using SpectralFitting.jl framework. Kept around so we can borrow the model data downloading functions.\n\n\n\n\n\n","category":"type"},{"location":"private/#BayesJ.call_mekal-Tuple{Any, Any, Any}","page":"Internals","title":"BayesJ.call_mekal","text":"call_mekal(energy_range, temperature, nH)\n\nGiven a unitful range of energy and unitless temperature (keV) and hydrogen density in the source (cm^-3) calls MEKAL to calculate the volume emissivity of the source in the source frame.\n\nDue to performance concerns I recommend calling prepare_model_mekal instead to generate an interpolation over the model. The interpolation also implicitly includes absorption, which this function does not.\n\nImplementation Details\n\nThe XSPEC and MEKAL documentation is at times short on detail. Here's what I've figured out, mostly by inspecting the models in the XSPEC source code. The MEKA model has similar code and more extensive comments that MEKAL  that are very helpful.\n\nThe values calculated by MEKAL are multiplied by a constant of 2.53325e-3. Inspecting the MEKA comments we find the following equation for this constant which we will term C\n\nC = 2  h^2  α^3  ((3πm_e)^15  (1000e)  π  (1pc)^2  10^4)\n\nand the instruction to use SI units. h appears to be Planck's constant (unreduced) and α the fine-structure constant. Calculating this gives C=2.533242122484874e-59 kg^1/2 m^2 A^-1/2 s^-5/2.\n\nMEKAL adjusts emissions by what it terms cem which is the emitting volume divided by the square of distance to the volume. This follows naturally from the relationship between surface emissivity σ and volume emissivity ϵ.\n\n    σ = ϵV  (4πD)^2\n\nMEKAL uses units of 1e50cm^3 / 1pc^2 for cem. If we multiply by 10^50cm^3 we get 2.533242122484874e-13 J^2 cm s^2 C^-1/2 kg^-3/2. As dividing by 1pc^2 would change more than just the magnitude of the value I assume it is folded in as the 1pc^2 in the original equation. Furthermore I assume that remaining difference rests in the details of units and implementation.\n\nMEKAL states it outputs in phot/cm^2/s/keV. MEKA states the same but has an optional switch that instead outputs in phot/m^3/s/keV. It does this by replacing C with D=3.03103E-9 which it states is derived so that\n\nC = D * 10^40  (4π * 1pc^2)\n\nI have verified this gives results matching that of doing the conversion manually (by dividing out 1e50cm^3 and multiplying by 1pc^2).\n\nAttempts were made to invoke MEKAL through XSPEC's wrapper functions but I found the modifications they introduced were troublesome to work out, though some were necessary and are reimplemented in surface_brightness. MEKAL's operations are more clearly physically motivated.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.prepare_model_mekal-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{Union{Unitful.Quantity{T, 𝐋^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-2, U}} where {L, S}} where {T, U}, AbstractRange{T}, Real}} where {T<:(Union{Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}} where {L, S}} where {T, U}), U<:(Union{Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}} where {L, S}} where {T, U}), V<:(Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U})}","page":"Internals","title":"BayesJ.prepare_model_mekal","text":"prepare_model_mekal(nHcol, energy_bins, z; temperatures, hydrogen_densities, use_interpolation=true)\n\nCreate an interpolated alias to the mekal model with photoelectric absorption over specified parameter ranges.\n\nThe interpolation object (T::Unitful.Energy, nH::Unitful.𝐋^-3) -> volume emissivity returned gives the volume emissivity per bin of a source of the specified temperature and hydrogen number density with the effects of absorption from passage through an area with a hydrogen column density equal to that specified in the prepare_model_mekal call.\n\nApplies redshift and time dilation.\n\nInterpolation has a significant performance improvement over calling the model directly. No redshift is currently applied to energy bins - they should be assumed to be in the source frame. It may be wise to adjust this but then surface_brightness will need adjustment so it doesn't apply the correction twice.\n\n\n\n\n\n","category":"method"},{"location":"private/#IO","page":"Internals","title":"IO","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"io.jl\", \"fitsio_fix.jl\"]","category":"page"},{"location":"private/#BayesJ.Dataset","page":"Internals","title":"BayesJ.Dataset","text":"Abstract type to wrap various specific formats for source data. By doing so we can make use of     multiple dispatch to automatically select the appropriate loading functions for data.\n\n\n\n\n\n","category":"type"},{"location":"private/#BayesJ.Ellipse","page":"Internals","title":"BayesJ.Ellipse","text":"Ellipse(x, y, r1, r2, θ)\n\nAn ellipse centered at (x, y) with radii r1 and r2 and rotation θ.\n\nr1 and r2 are the width and height respectively with zero rotation and θ is defined as rotation counterclockwise from the x-axis in radians.\n\n\n\n\n\n","category":"type"},{"location":"private/#BayesJ.load_data-Tuple{BayesJ.Dataset}","page":"Internals","title":"BayesJ.load_data","text":"load_data(data)\n\nLoad events data from a given dataset, returning pairs of events and exposure times.\n\nThe first pair includes all observed events, the second is the background observation.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.load_mask-Tuple{AbstractString, Any, Any}","page":"Internals","title":"BayesJ.load_mask","text":"load_mask(path, x_edges, y_edges)\n\nTakes a mask file as used by DS9 and x and y edges to bin on. This should match the edges used to bin the data.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.load_mask-Tuple{BayesJ.Dataset}","page":"Internals","title":"BayesJ.load_mask","text":"load_mask(data)\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.load_response-Tuple{BayesJ.Dataset, Any}","page":"Internals","title":"BayesJ.load_response","text":"load_response(data, energy_range)\n\nLoad the RMF and ARF for an observation, trimmed for the desired energy range\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.test_point","page":"Internals","title":"BayesJ.test_point","text":"test_point(e::Ellipse, x, y, atol=0)\n\nTests whether a point is within an ellipse with tolerance atol.\n\nThe equation for an ellipse in Cartesian coordinates is of the form f(x,y) = 1. We consider any point (x0, y0) such that f(x0, y0) <= 1 + atol as being within the ellipse.\n\n\n\n\n\n","category":"function"},{"location":"private/#Likelihood","page":"Internals","title":"Likelihood","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"likelihood.jl\"]","category":"page"},{"location":"private/#BayesJ.Prior","page":"Internals","title":"BayesJ.Prior","text":"Abstract supertype for priors. Should implement:\n\nA transform(prior, x) function that transforms a value x on the unit range to a value on the distribution represented by the prior.\nA field name<:AbstractString \n\n\n\n\n\n","category":"type"},{"location":"private/#BayesJ.log_factorial-Tuple{N} where N<:Integer","page":"Internals","title":"BayesJ.log_factorial","text":"log_factorial(n)\n\nFinds lnn the natural logarithm of the factorial of n.\n\nn rapidly gets to large to quickly and directly calculate the factorial so we exploit logarithm rules to expand it out to a series of sums.\n\nIt is intended to be broadcast across all values of the data array.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.log_likelihood-NTuple{5, Any}","page":"Internals","title":"BayesJ.log_likelihood","text":"log_likelihood(observed, observed_background, predicted, predicted_background, observed_log_factorial)\n\nCalculate the log-likelihood of the prediction given an observation.\n\nThe observed and predicted arrays include background events. Log factorial is calculated as ln(observed) + ln(observed_background). We require it to be supplied to improve performance - no need to calculate it every time.\n\nWe assume the predicted_background is scaled to the same exposure time as the observed background.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.make_cube_transform-Tuple{Vararg{BayesJ.Prior}}","page":"Internals","title":"BayesJ.make_cube_transform","text":"make_cube_transform(priors::Prior...)::NTuple{2, Function}\n\nTurn a sequence of prior transform functions into a transform function transform_cube that operates on the hypercube generated by multinest, stripping delta functions.\n\nIt also returns a function reconstruct_args take takes the output of the first function and returns it with the values of the delta priors inserted in appropriate postions.  This is because Ultranest can run with delta priors but outputs error warnings in the process.\n\nReturns transform_cube and reconstruct_args as a tuple.\n\n\n\n\n\n","category":"method"},{"location":"private/#BayesJ.transform-Tuple{BayesJ.Prior, Real}","page":"Internals","title":"BayesJ.transform","text":"transform(prior, x)\n\nTransforms a value x on the unit range to a value on the distribution represented by the prior.\n\n\n\n\n\n","category":"method"},{"location":"private/#Utilities","page":"Internals","title":"Utilities","text":"","category":"section"},{"location":"private/","page":"Internals","title":"Internals","text":"Modules = [BayesJ]\nPublic  = false\nPrivate = true\nPages   = [\"mpi.jl\", \"run.jl\"]","category":"page"},{"location":"public/#Public-Documentation","page":"Public","title":"Public Documentation","text":"","category":"section"},{"location":"public/#Index","page":"Public","title":"Index","text":"","category":"section"},{"location":"public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Public-Interface","page":"Public","title":"Public Interface","text":"","category":"section"},{"location":"public/","page":"Public","title":"Public","text":"Modules = [BayesJ]\nPublic  = true\nPrivate = false","category":"page"},{"location":"public/#BayesJ.DeltaPrior","page":"Public","title":"BayesJ.DeltaPrior","text":"DeltaPrior(name::AbstractString, value::Number)\n\nA delta prior that always returns a constant value.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.FITSData","page":"Public","title":"BayesJ.FITSData","text":"FITSData(observation::AbstractString, background::AbstractString, arf::AbstractString, rmf::AbstractString, pixel_edge_angle::DimensionfulAngles.Angle)\n\nStore a collection of observation data from fits files.\n\nTakes paths to the files as inputs.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.GenericPrior","page":"Public","title":"BayesJ.GenericPrior","text":"GenericPrior(name::AbstractString, dist::Distributions.UnivariateDistribution)\n\nA generic prior that draws from the specified distribution.\n\nSee Distributions.jl for a list of avaliable distributions.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.LogUniformPrior","page":"Public","title":"BayesJ.LogUniformPrior","text":"LogUniformPrior(name::AbstractString, min::Number, max::Number)\n\nA log uniform prior that draws from a distribution between min and max whose base 10 logarithm is uniformly distributed.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.NormalPrior","page":"Public","title":"BayesJ.NormalPrior","text":"NormalPrior(name::AbstractString, mean::Number, σ::Number)\n\nA normal prior that draws from a Gaussian/normal distribution with mean and standard deviation σ.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.PlaintextData","page":"Public","title":"BayesJ.PlaintextData","text":"PlaintextData(observation::AbstractString, background::AbstractString, arf::AbstractString, rmf::AbstractString, shape::NTuple{3,Int}, rmf_shape::NTuple{2,Int}, obs_exposure_time::Unitful.Time, bg_exposure_time::Unitful.Time, pixel_edge_angle::DimensionfulAngles.Angle, rmf_unit::Unitful.Area)\n\nStore preprocessed data from plaintext files with whitespace delimiters.\n\nEvents and response files provided as paths. This data is assumed to have the appropriate binning and energy ranges already applied. The shape of the event arrays should be (channels, x, y) and the rmf_shape should be (channels, bins).\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.UniformPrior","page":"Public","title":"BayesJ.UniformPrior","text":"UniformPrior(name::AbstractString, min::Number, max::Number)\n\nA uniform prior that draws from a uniform distribution between min and max.\n\n\n\n\n\n","category":"type"},{"location":"public/#BayesJ.Model_Einasto-Tuple{Real, Vararg{Any, 6}}","page":"Public","title":"BayesJ.Model_Einasto","text":"Mass is in solar masses.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_Einasto-Tuple{Union{Unitful.Quantity{T, 𝐌, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐌, U}} where {L, S}} where {T, U}, Vararg{Any, 6}}","page":"Public","title":"BayesJ.Model_Einasto","text":"Model_Einasto(MT_200, fg_200, α, a, b, c, c_500_GNFW, z)\n\nGenerate a cluster profile based on the Einasto mass density profile and the GNFW gas pressure profile.\n\nReturns functions for gas temperature and gas mass density as a function of radius. See Appendix A, Olamaie et al. 2015 (doi:10.1093/mnras/stu2146)\n\nIf α>2 then the gas temperature starts increasing to physically improbable levels (orders exceeding 10¹¹ keV). We thus constrain the value to below that. \n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_NFW-Tuple{Real, Vararg{Any, 5}}","page":"Public","title":"BayesJ.Model_NFW","text":"Model_NFW(MT_200::Real, fg_200, a, b, c, c_500_GNFW, z)\n\nGenerate a cluster profile based on the NFW mass density and GNFW gas density profiles.\n\nMass is in solar masses.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_NFW-Tuple{Union{Unitful.Quantity{T, 𝐌, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐌, U}} where {L, S}} where {T, U}, Vararg{Any, 5}}","page":"Public","title":"BayesJ.Model_NFW","text":"Model_NFW(MT_200::Unitful.Mass, fg_200, a, b, c, c_500_GNFW, z)\n\nGenerate a cluster profile based on the NFW mass density and GNFW gas density profiles.\n\nUses the model from Olamaie 2012 (doi:10.1111/j.1365-2966.2012.20980.x), which is based on the NFW dark matter density profile and the GNFW gas pressure profile. Returns functions for gas temperature and gas mass density as a function of radius.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.Model_Vikhlinin2006-Tuple{Real, Real, Real, Real, Any, Any, Any, Any, Real, Real, Any, Real, Any, Real, Any, Any, Any}","page":"Public","title":"BayesJ.Model_Vikhlinin2006","text":"Densities are in cm^-3 and radii are in kpc. Temperature in keV.\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.sample-Tuple{BayesJ.Dataset, AbstractRange{<:Union{Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^2 𝐌 𝐓^-2, U}} where {L, S}} where {T, U}}, Function, AbstractVector{<:BayesJ.Prior}, Union{Unitful.Quantity{T, 𝐋^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-2, U}} where {L, S}} where {T, U}, Real, Tuple{var\"#s1796\", var\"#s1796\"} where var\"#s1796\"<:Real, Tuple{var\"#s1797\", var\"#s1797\"} where var\"#s1797\"<:Real}","page":"Public","title":"BayesJ.sample","text":"sample(data::Dataset, energy_range::AbstractRange{Unitful.Energy}, cluster_model::Function, priors::AbstractVector{Prior}, nhCol::SurfaceDensity, redshift, x, y)\n\nRun Bayesian inference on a given set of data considering only the selected energy range.\n\nAn gas emission model (density, temperature) → emissivity can be provided.\nThe first two priors should always be x0 and y0, giving cluster centre position.\nx and y are tuples of (min, max).\n\n\n\n\n\n","category":"method"},{"location":"public/#BayesJ.sample-Union{Tuple{T}, Tuple{T, T, Matrix, Function, Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐓, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐓, U}} where {L, S}} where {T, U}, Real}} where T<:AbstractArray","page":"Public","title":"BayesJ.sample","text":"sample(observed, observed_background, response_function, transform, obs_exposure_time, bg_exposure_time, redshift; prior_names, cluster_model, emission_model, param_wrapper, pixel_edge_angle, background_rate, average_effective_area)\n\nConfigure some necessary variables and launch ultranest.\n\nThe observed array includes the background.\nThe response function includes both the RMF and ARF, as described in apply_response_function.\nThe emission model should be a function compatible with the requirements of the surface_brightness function, which it will be passed to.\nThe pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.\n\nThis area is assumed to be square with the edge angle giving the side length.\n\nThe average effective area is the effective area of the telescope averaged across energies,\n\nused with the total background rate across all channels (counts per unit telescope area per sky angle per second)  to calculate the background counts per second per channel per pixel.\n\nThe first two priors should always be \"x0\" and \"y0\", giving cluster centre position and the order of prior_names must match the transform.\nThe cluster model should be a function that takes the parameters (and redshift as a kwarg) and returns (gas_temperature,gas_density) as functions of\n\nradius which return their respective quantities with units.\n\nparam_wrapper takes the output of the transform function and adds any additional arguements necessary for the model.\n\n\n\n\n\n","category":"method"},{"location":"#BayesJ.jl-Documentation","page":"Home","title":"BayesJ.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
