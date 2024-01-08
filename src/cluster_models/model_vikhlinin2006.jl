export Model_Vikhlinin2006

"""
    Model_Vikhlinin2006(n0::NumberDensity, n02::NumberDensity, rc::Unitful.Length, rc2::Unitful.Length, α, β, β2, ϵ, rs::Unitful.Length, T0::Unitful.Energy, TminT0, rcool::Unitful.Length, acool, rt::Unitful.Length, a, b, c; γ=3, kwargs...)

Generate a cluster profile using the highly free models from [vikhlininChandraSampleNearby2006](@cite).


"""
function Model_Vikhlinin2006(
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
    c,
    d;
    γ=3,
    kwargs...
)::NTuple{2,Function}

    priorcheck(ϵ < 5, -1e100(1 + ϵ))

    # Performance of exponents sucks if we don't strip units
    n0_u::Float64 = ustrip(Float64, u"cm^-3", n0)
    n02_u::Float64 = ustrip(Float64, u"cm^-3", n02)
    rc_u::Float64 = ustrip(Float64, u"kpc", rc)
    rc2_u::Float64 = ustrip(Float64, u"kpc", rc2)
    rs_u::Float64 = ustrip(Float64, u"kpc", rs)

    priorcheck(rc2_u < rc_u, -1e100 * (1 + (rc2_u - rc_u)))
    priorcheck(rc_u < rs_u, -1e100 * (1 + (rc_u - rs_u)))

    @mpidebug "Model Vikhlinin2006 called with" n0_u n02_u rc_u rc2_u α β β2 ϵ rs_u T0 TminT0 rcool acool rt a b c γ

    """
        np_ne(r, n0, n02, rc, rc2, α, β, β2, ϵ, rs; γ=3)

    Equation 3 from Vikhlinin et al. 2006. Calculates emission measure.

    They constain ϵ<5 to exclude unphysically sharp density breaks.
    """
    function np_ne(
        r,
        n0,
        n02,
        rc,
        rc2,
        α,
        β,
        β2,
        ϵ,
        rs;
        γ=3
    )
        n0^2 *
        (r / rc)^(-α) / (1 + r^2 / rc^2)^(3β - α / 2) *
        1 / (1 + r^γ / rs^γ)^(ϵ / γ) +
        n02^2 / (1 + r^2 / rc2^2)^(3β2)
    end

    """
        gas_density(r, ne_np)

    Calculates gas density from the emission measure.
    """
    function gas_density(r::Unitful.Length)::Unitful.Density
        r = abs(r)
        res = let n0 = n0_u, n02 = n02_u, rc = rc_u, rc2 = rc2_u, α = α, β = β, β2 = β2, ϵ = ϵ, rs = rs_u, γ = γ
            μ * 1u"cm^-3" * sqrt(
                np_ne(ustrip(u"kpc", r), n0, n02, rc, rc2, α, β, β2, ϵ, rs, γ=γ)
            )
        end
        # We round to zero because it keeps MEKAL behaving
        isapprox(res, 0u"g/cm^3", atol=1e-30u"g/cm^3") ? 0.0u"g/cm^3" : res
    end

    T0 = ustrip(u"keV", T0)
    rcool = ustrip(u"kpc", rcool)
    rt = ustrip(u"kpc", rt)

    """
        t(r, rt, a, b, c)

    Equation 4 from Vikhlinin et al. 2006. Models temperature profile outside
    central cooling region.
    """
    function t(r, rt, a, b, c, d)
        # Tweaking Vikhlinin's equation to avoid an infinite temperature at the core
        (r / rt + d)^(-a) / (1 + (r / rt)^b)^(c / b)
    end

    """
        tcool(r, rcool, acool, Tmin, T0)

    Equation 5 from Vikhlinin et al. 2006. Models temperature profile in 
    cluster core.

    TminT0 = Tmin/T0
    """
    function tcool(r, rcool, acool, TminT0)
        x = (r / rcool)^acool
        (x + TminT0) / (x + 1)
    end

    """
        T3D(r, T0, Tmin, rcool, acool, rt, a, b, c)

    Equation 6 from Vikhlinin et al. 2006. Combines [`t`](@ref) and [`tcool`](@ref)
    to model the temperature profile throughout the cluster.
    """
    function gas_temperature(r, T0, TminT0, rcool, acool, rt, a, b, c, d)
        T0 *
        tcool(r, rcool, acool, TminT0) *
        t(r, rt, a, b, c, d)
    end
    function gas_temperature(r::Unitful.Length)::Unitful.Energy
        r = abs(r)
        res = let T0 = T0, TminT0 = TminT0, rcool = rcool, acool = acool, rt = rt, a = a, b = b, c = c, d = d
            1u"keV" * gas_temperature(
                ustrip(u"kpc", r), T0, TminT0, rcool, acool, rt, a, b, c, d
            )
        end
        isapprox(res, 0u"keV", atol=1e-30u"keV") ? 0.0u"keV" : res
    end


    # if !isfinite(gas_temperature(0u"kpc"))
    #     @mpirankedwarn "Infinite temperature at core" rt a b c (0 / rt) (0 / rt)^(-a)
    #     throw(PriorError(-2e100))
    # end

    # @assert isfinite(gas_density(1u"kpc")) "Gas density not finite: $([n0_u, n02_u, rc_u, rc2_u, α, β, β2, ϵ, rs_u, γ])"
    # @assert isfinite(gas_temperature(1u"kpc")) "Gas temperature not finite: $([T0, TminT0, rcool, acool, rt, a, b, c])"

    return gas_temperature, gas_density

end
"""
    Model_Vikhlinin2006(n0::Real, n02::Real, rc::Real, rc2::Real, α, β, β2, ϵ, rs::Real, T0::Real, TminT0, rcool::Real, acool, rt::Real, a, b, c; γ=3, kwargs...)

Unitless wrapper for [`Model_Vikhlinin2006`](@ref). Densities are in cm^-3 and radii are in kpc. Temperature is in keV.
"""
function Model_Vikhlinin2006(
    n0::Real,
    n02::Real,
    rc::Real,
    rc2::Real,
    α,
    β,
    β2,
    ϵ,
    rs::Real,
    T0::Real,
    TminT0,
    rcool::Real,
    acool,
    rt::Real,
    a,
    b,
    c,
    d;
    γ=3,
    kwargs...
)::NTuple{2,Function}
    Model_Vikhlinin2006(
        n0 * 1u"cm^-3",
        n02 * 1u"cm^-3",
        rc * 1u"kpc",
        rc2 * 1u"kpc",
        α,
        β,
        β2,
        ϵ,
        rs * 1u"kpc",
        T0 * 1u"keV",
        TminT0,
        rcool * 1u"kpc",
        acool,
        rt * 1u"kpc",
        a,
        b,
        c,
        d;
        γ=γ
    )
end
