using Cosmology
using Unitful
using PhysicalConstants.CODATA2018: m_p

const cosmo::Cosmology = cosmology()
const μ = 0.6 * m_p
const μ_e = 1.14 * m_p