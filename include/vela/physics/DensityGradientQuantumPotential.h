#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/mesh/DeviceMesh.h"

#include <unordered_map>
#include <string>
#include <vector>

namespace vela {

/// Sentaurus-compatible electron density-gradient correction for Boltzmann
/// statistics.  The solved energy divided by q is stored in volts.
struct DensityGradientQuantumPotentialConfig {
    bool enabled = false;
    std::string couplingMode = "outer"; ///< "outer" or frozen imported potential.
    Real gamma = 3.6;
    Real effectiveMassRatio = 1.0618016171622988;
    int maxIterations = 30;
    Real relativeTolerance = 1.0e-7;
    Real absoluteTolerance_V = 1.0e-10;
    Real damping = 0.5;
    Real maxUpdate_V = 0.1;
    int outerMaxIterations = 20;
};

struct DensityGradientQuantumPotentialResult {
    VectorXd potential_V;
    int iterations = 0;
    Real residualInfinityNorm = 0.0;
    bool converged = false;
};

/// Return gamma*hbar^2/(6*m*q), in V*m^2.
Real densityGradientCoefficientVm2(
    const DensityGradientQuantumPotentialConfig& config);

/// Solve Lambda*u + C*laplacian(u)=0 with
/// u=sqrt(n/Nref) and n=classicalDensity*exp(-Lambda/Vt).
/// Edges to inactive (insulating) nodes are omitted, giving the homogeneous
/// Neumann boundary used at semiconductor/insulator interfaces.  Values on
/// nodes listed in dirichletPotential_V are imposed exactly.
DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units = PhysicalUnitSystem::legacySI(),
    DensityGradientQuantumPotentialConfig config = {},
    const VectorXd& initialPotential_V = {});

} // namespace vela
