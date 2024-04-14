export Vikhlinin2006Model

"""
    Vikhlinin2006Model(n0::N, n02::N, rc::L, rc2::L, α, β, β2, ϵ, rs::L, T0::T, TminT0, rcool::L, acool, rt::L, a, b, c, γ=3) where {N<:NumberDensity,L<:Unitful.Length,T<:Unitful.Energy}

Generate a cluster profile using the highly free models from [vikhlininChandraSampleNearby2006](@cite).
"""
struct Vikhlinin2006Model{N<:NumberDensity,L<:Unitful.Length,T<:Unitful.Energy} <: ClusterModel
    n0::N
    n02::N
    rc::L
    rc2::L
    α
    β
    β2
    ϵ
    rs::L
    T0::T
    TminT0
    rcool::L
    acool
    rt::L
    a
    b
    c
    γ
end


"""
    Vikhlinin2006Model(n0::NumberDensity, n02::NumberDensity, rc::Unitful.Length, rc2::Unitful.Length, α, β, β2, ϵ, rs::Unitful.Length, T0::Unitful.Energy, TminT0, rcool::Unitful.Length, acool, rt::Unitful.Length, a, b, c; γ=3, kwargs...)


Constructor for [`Vikhlinin2006Model`](@ref) that performs prior checks.

If not specified densities are assumed to be in cm^-3, radii are in kpc and temperature is in keV.
"""
function Vikhlinin2006Model(
    n0::NumberDensity,
    n02::NumberDensity,
    rc::Unitful.Length,
    rc2::Unitful.Length,
    α,
    β,
    β2,
    ϵ,
    rs::Unitful.Length,
    T0::Unitful.Energy,
    TminT0,
    rcool::Unitful.Length,
    acool,
    rt::Unitful.Length,
    a,
    b,
    c;
    γ=3,
    kwargs...
)

    priorcheck(ϵ < 5, -1e100(1 + ϵ))

    # Performance of exponents sucks if we don't strip units
    # n0_u::Float64 = ustrip(Float64, u"cm^-3", n0)
    # n02_u::Float64 = ustrip(Float64, u"cm^-3", n02)
    # rc_u::Float64 = ustrip(Float64, u"kpc", rc)
    # rc2_u::Float64 = ustrip(Float64, u"kpc", rc2)
    # rs_u::Float64 = ustrip(Float64, u"kpc", rs)

    priorcheck(rc2 < rc, -1e100 * (1 + ustrip(u"Mpc", rc2 - rc)))
    priorcheck(rc < rs, -1e100 * (1 + ustrip(u"Mpc", rc - rs)))

    # @mpidebug "Model Vikhlinin2006 called with" n0_u n02_u rc_u rc2_u α β β2 ϵ rs_u T0 TminT0 rcool acool rt a b c γ
    return Vikhlinin2006Model(n0, n02, uconvert(u"kpc", rc), uconvert(u"kpc", rc2), α, β, β2, ϵ, uconvert(u"kpc", rs), T0, TminT0, uconvert(u"kpc", rcool), acool, uconvert(u"kpc", rt), a, b, c, γ)
end
Vikhlinin2006Model(n0::Real, n02::Real, rc::Real, rc2::Real, α, β, β2, ϵ, rs::Real, T0::Real, TminT0, rcool::Real, acool, rt::Real, a, b, c; γ=3, kwargs...) = Vikhlinin2006Model(n0 * 1u"cm^-3", n02 * 1u"cm^-3", rc * 1u"kpc", rc2 * 1u"kpc", α, β, β2, ϵ, rs * 1u"kpc", T0 * 1u"keV", TminT0, rcool * 1u"kpc", acool, rt * 1u"kpc", a, b, c; γ=γ)

"""
    np_ne(cm::Vikhlinin2006Model)

Equation 3 from Vikhlinin et al. 2006. Calculates emission measure.

They constain ϵ<5 to exclude unphysically sharp density breaks.
"""
function np_ne(cm::Vikhlinin2006Model, r::Unitful.Length{<:Real})
    cm.n0^2 *
    (r / cm.rc)^(-cm.α) / (1 + r^2 / cm.rc^2)^(3cm.β - cm.α / 2) *
    1 / (1 + r^cm.γ / cm.rs^cm.γ)^(cm.ϵ / cm.γ) +
    cm.n02^2 / (1 + r^2 / cm.rc2^2)^(3cm.β2)
end

"""
    density(cm::Vikhlinin2006Model, r::Unitful.Length{<:Real})

Calculates gas density from the emission measure.
"""
function density(cm::Vikhlinin2006Model, r::Unitful.Length{<:Real})
    r = abs(r)
    res = μ * sqrt(np_ne(cm, r))
    # We round to zero because it keeps MEKAL behaving
    isapprox(res, 0u"g/cm^3", atol=1e-30u"g/cm^3") ? 0.0u"g/cm^3" : res
end

"""
    temperature(cm::Vikhlinin2006Model, r::Unitful.Length{<:Real})

Equation 6 from Vikhlinin et al. 2006.
"""
function temperature(cm::Vikhlinin2006Model, r::Unitful.Length{<:Real})
    r = abs(r)
    x = (r / cm.rcool)^cm.acool
    res = cm.T0 * (x + cm.TminT0) / (x + 1) * (r / cm.rt)^(-cm.a) / (1 + (r / cm.rt)^cm.b)^(cm.c / cm.b)
    isapprox(res, 0u"keV", atol=1e-30u"keV") ? 0.0u"keV" : res
end