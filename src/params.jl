using Cosmology
using Unitful
using PhysicalConstants.CODATA2018: m_p

"""Cosmology assumed by BayesX"""
const cosmo = cosmology() # todo: make configurable

"""Mean mass per gas particle"""
const μ = 0.6 * m_p

"""Mean gas mass per electron"""
const μ_e = 1.14 * m_p