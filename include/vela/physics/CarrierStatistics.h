#pragma once

#include "vela/material/Material.h"
#include <string>

namespace vela {

struct CarrierStatisticsConfig {
    /// "boltzmann" or "fermi_dirac".  Boltzmann remains the compatibility default.
    std::string model = "boltzmann";
};

struct EquilibriumCarrierState {
    Real potential = 0.0; ///< Electrostatic potential relative to the intrinsic reference [V].
    Real n = 0.0;         ///< Electron density [m^-3].
    Real p = 0.0;         ///< Hole density [m^-3].
};

bool usesFermiDirac(const CarrierStatisticsConfig& config);

/**
 * @brief Normalized complete Fermi-Dirac integral F_{1/2}(eta).
 *
 * Uses the Bednarczyk analytic approximation (maximum relative error 0.4%).
 * The normalization is 2/sqrt(pi) times the defining integral, so the
 * non-degenerate limit is exp(eta).
 */
Real fermiDiracHalf(Real eta);
Real fermiDiracHalfDerivative(Real eta);
Real inverseFermiDiracHalf(Real value);

/**
 * @brief Boltzmann carrier statistics.
 *
 * Boltzmann (non-degenerate) approximation:
 *   n = ni * exp((psi - phin) / Vt)
 *   p = ni * exp((phip - psi) / Vt)
 *
 * Exponent arguments are clamped to [-500, 500] to prevent overflow.
 *
 * @param ni    Intrinsic carrier concentration [m^-3]
 * @param psi   Electrostatic potential [V]
 * @param phin  Electron quasi-Fermi potential [V]
 * @param phip  Hole quasi-Fermi potential [V]
 * @param Vt    Thermal voltage kT/q [V]
 */
double electronDensity(double ni, double psi, double phin, double Vt);
double holeDensity    (double ni, double psi, double phip, double Vt);

/// Carrier densities with an explicit statistics model and density of states.
Real electronDensity(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                     const CarrierStatisticsConfig& config);
Real holeDensity(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                 const CarrierStatisticsConfig& config);

/// Derivatives with respect to the dimensionless reduced Fermi energy.
Real electronDensityDerivativeEta(Real ni, Real Nc, Real psi, Real phin, Real Vt,
                                  const CarrierStatisticsConfig& config);
Real holeDensityDerivativeEta(Real ni, Real Nv, Real psi, Real phip, Real Vt,
                              const CarrierStatisticsConfig& config);

Real electronQuasiFermiPotential(Real ni, Real Nc, Real psi, Real n, Real Vt,
                                 const CarrierStatisticsConfig& config);
Real holeQuasiFermiPotential(Real ni, Real Nv, Real psi, Real p, Real Vt,
                             const CarrierStatisticsConfig& config);

/**
 * @brief Charge-neutral equilibrium state used by ideal Ohmic contacts.
 *
 * Solves n(psi)-p(psi)=netDoping numerically for Fermi-Dirac statistics.
 * The Boltzmann branch retains the closed-form legacy solution.
 */
EquilibriumCarrierState equilibriumCarrierState(
    Real netDoping, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config);

/// Equilibrium n0*p0 at the same local charge imbalance n-p.
Real equilibriumCarrierProduct(
    Real n, Real p, Real ni, Real Nc, Real Nv, Real Vt,
    const CarrierStatisticsConfig& config);

/// Temperature-adjusted intrinsic density for a material using temperature_K.
double intrinsicDensity(const Material& material, double temperature_K);

} // namespace vela
