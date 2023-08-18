export Model_Vikhlinin2006

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
    c;
    γ=3,
    kwargs...
)::NTuple{2,Function}

    # Performance of exponents sucks if we don't strip units
    n0 = ustrip(u"cm^-3", n0)
    n02 = ustrip(u"cm^-3", n02)
    rc = ustrip(u"kpc", rc)
    rc2 = ustrip(u"kpc", rc2)
    rs = ustrip(u"kpc", rs)

    @mpidebug "Model Vikhlinin2006 called with" n0 n02 rc rc2 α β β2 ϵ rs T0 TminT0 rcool acool rt a b c γ

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

    Vikhlinin et al. 2006 says
    > For the cosmic plasma with primordial He abundance and
    > abundances of heavier elements Z = 0.2 Z⊙
    > ``ρ_g = 1.624 m_p (n_p n_e)^{1/2}``
    I'm trying 1.14mₚ instead to match our other models
    and the metallicity assumptions in MEKAL.
    """
    function gas_density(r::Unitful.Length)::Unitful.Density
        r = abs(r)
        res = let
            n0 = n0
            n02 = n02
            rc = rc
            rc2 = rc2
            α = α
            β = β
            β2 = β2
            ϵ = ϵ
            rs = rs
            γ = γ
            μ_e * 1u"cm^-3" * sqrt(
                np_ne(ustrip(u"kpc", r), n0, n02, rc, rc2, α, β, β2, ϵ, rs, γ=γ)
            )
        end
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
    function t(r, rt, a, b, c)
        (r / rt)^(-a) / (1 + (r / rt)^b)^(c / b)
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
    function gas_temperature(r, T0, TminT0, rcool, acool, rt, a, b, c)
        T0 *
        tcool(r, rcool, acool, TminT0) *
        t(r, rt, a, b, c)
    end
    function gas_temperature(r::Unitful.Length)::Unitful.Energy
        r = abs(r)
        res = let
            T0 = T0
            TminT0 = TminT0
            rcool = rcool
            acool = acool
            rt = rt
            a = a
            b = b
            c = c
            1u"keV" * gas_temperature(
                ustrip(u"kpc", r), T0, TminT0, rcool, acool, rt, a, b, c
            )
        end
        isapprox(res, 0u"keV", atol=1e-30u"keV") ? 0.0u"keV" : res
    end

    @assert isfinite(gas_density(1u"kpc")) "Gas density not finite: $([n0, n02, rc, rc2, α, β, β2, ϵ, rs, γ])"
    @assert isfinite(gas_temperature(1u"kpc")) "Gas temperature not finite: $([T0, TminT0, rcool, acool, rt, a, b, c])"

    return gas_temperature, gas_density

end
"""
Densities are in cm^-3 and radii are in kpc. Temperature in keV.
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
    c;
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
        c;
        γ=γ
    )
end
