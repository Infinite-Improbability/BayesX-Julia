using Cosmology
using Unitful
using PhysicalConstants.CODATA2018: m_p, G, k_B

"""
Cosmology assumed by BayesX

Defaults to:
* h=0.69,
* Neff=3.04,
* OmegaK=0,
* OmegaM=0.29,
* OmegaR=nothing,
* Tcmb=2.7255,
* w0=-1,
* wa=0
"""
const cosmo = cosmology(h=0.69,
    Neff=3.04,
    OmegaK=0,
    OmegaM=0.29,
    OmegaR=nothing,
    Tcmb=2.7255,
    w0=-1,
    wa=0) # todo: make configurable

"""Mean mass per gas particle = 0.6mₚ"""
const μ = 0.6 * m_p

"""Mean gas mass per electron = 1.14mₚ"""
const μ_e = 1.14 * m_p