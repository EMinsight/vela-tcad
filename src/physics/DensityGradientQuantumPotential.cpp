#include "vela/physics/DensityGradientQuantumPotential.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/solver/LinearSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
namespace {

Real safeExponential(Real argument)
{
    constexpr Real limit = 700.0;
    return std::exp(std::clamp(argument, -limit, limit));
}

} // namespace

Real densityGradientCoefficientVm2(
    const DensityGradientQuantumPotentialConfig& config)
{
    if (!(config.gamma > 0.0) || !std::isfinite(config.gamma) ||
        !(config.effectiveMassRatio > 0.0) ||
        !std::isfinite(config.effectiveMassRatio)) {
        throw std::invalid_argument(
            "density-gradient gamma and effective mass must be finite and positive.");
    }
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    return config.gamma * hbar * hbar /
        (6.0 * config.effectiveMassRatio * constants::m0 * constants::q);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    const int nodeCount = static_cast<int>(mesh.numNodes());
    if (classicalDensity.size() != nodeCount ||
        activeNodes.size() != static_cast<std::size_t>(nodeCount)) {
        throw std::invalid_argument(
            "density-gradient input sizes must match the mesh node count.");
    }
    if (!(thermalVoltage_V > 0.0) || !std::isfinite(thermalVoltage_V))
        throw std::invalid_argument("density-gradient thermal voltage must be positive.");
    if (config.maxIterations <= 0 || !(config.relativeTolerance >= 0.0) ||
        !(config.absoluteTolerance_V >= 0.0) || !(config.damping > 0.0) ||
        config.damping > 1.0 || !std::isfinite(config.damping) ||
        !(config.maxUpdate_V >= 0.0) || !std::isfinite(config.maxUpdate_V)) {
        throw std::invalid_argument("invalid density-gradient iteration configuration.");
    }

    VectorXd potential = initialPotential_V.size() == nodeCount
        ? initialPotential_V : VectorXd::Zero(nodeCount);
    for (int i = 0; i < nodeCount; ++i) {
        if (!std::isfinite(classicalDensity(i)) || classicalDensity(i) < 0.0)
            throw std::invalid_argument("density-gradient density must be finite and non-negative.");
        if (!std::isfinite(potential(i)))
            potential(i) = 0.0;
        if (!activeNodes[static_cast<std::size_t>(i)])
            potential(i) = 0.0;
    }
    for (const auto& [node, value] : dirichletPotential_V) {
        if (node >= mesh.numNodes() || !std::isfinite(value))
            throw std::invalid_argument("invalid density-gradient Dirichlet boundary.");
        potential(static_cast<int>(node)) = value;
    }

    const Real coefficient = densityGradientCoefficientVm2(config);
    const Real lengthScale_m = units.lengthMPerInternal();
    const Real areaScale_m2 = lengthScale_m * lengthScale_m;
    const Real referenceDensity = std::max<Real>(classicalDensity.maxCoeff(), 1.0);
    const Real densityFloor = std::numeric_limits<Real>::min();
    LinearSolver solver;
    DensityGradientQuantumPotentialResult result;

    for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
        VectorXd rootDensity(nodeCount);
        VectorXd derivative(nodeCount);
        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                rootDensity(i) = 0.0;
                derivative(i) = 0.0;
                continue;
            }
            const Real normalized = std::max(classicalDensity(i), densityFloor) /
                referenceDensity;
            rootDensity(i) = std::sqrt(normalized) *
                safeExponential(-potential(i) / (2.0 * thermalVoltage_V));
            derivative(i) = -rootDensity(i) / (2.0 * thermalVoltage_V);
        }

        VectorXd residual = VectorXd::Zero(nodeCount);
        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(static_cast<std::size_t>(nodeCount) * 8);
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            if (!activeNodes[edge.n0] || !activeNodes[edge.n1] ||
                edge.length <= 1.0e-30 || edge.couple <= 0.0)
                continue;
            const Real laplaceWeight = coefficient * edge.couple / edge.length;
            residual(i) += laplaceWeight * (rootDensity(j) - rootDensity(i));
            residual(j) += laplaceWeight * (rootDensity(i) - rootDensity(j));
            triplets.emplace_back(i, i, -laplaceWeight * derivative(i));
            triplets.emplace_back(i, j,  laplaceWeight * derivative(j));
            triplets.emplace_back(j, i,  laplaceWeight * derivative(i));
            triplets.emplace_back(j, j, -laplaceWeight * derivative(j));
        }

        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                residual(i) = potential(i);
                triplets.emplace_back(i, i, 1.0);
                continue;
            }
            const Real volume_m2 = mesh.getNode(static_cast<Index>(i)).volume * areaScale_m2;
            residual(i) += potential(i) * rootDensity(i) * volume_m2;
            triplets.emplace_back(
                i, i, (rootDensity(i) + potential(i) * derivative(i)) * volume_m2);
        }

        for (const auto& [node, value] : dirichletPotential_V) {
            const int row = static_cast<int>(node);
            residual(row) = potential(row) - value;
            triplets.emplace_back(row, row, 1.0);
        }

        SparseMatrixd jacobian(nodeCount, nodeCount);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
        // Replace constrained rows after summation, including edge entries.
        for (const auto& [node, value] : dirichletPotential_V) {
            (void)value;
            const int row = static_cast<int>(node);
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        for (int row = 0; row < nodeCount; ++row) {
            if (activeNodes[static_cast<std::size_t>(row)])
                continue;
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        jacobian.prune(0.0);

        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        VectorXd update = solver.solve(jacobian, -residual);
        Real updateNorm = update.lpNorm<Eigen::Infinity>();
        if (config.maxUpdate_V > 0.0 && updateNorm > config.maxUpdate_V) {
            update *= config.maxUpdate_V / updateNorm;
            updateNorm = config.maxUpdate_V;
        }
        potential += config.damping * update;
        for (const auto& [node, value] : dirichletPotential_V)
            potential(static_cast<int>(node)) = value;
        result.iterations = iteration;
        const Real scale = std::max<Real>(potential.lpNorm<Eigen::Infinity>(), 1.0);
        if (updateNorm <= config.absoluteTolerance_V + config.relativeTolerance * scale) {
            result.converged = true;
            break;
        }
    }
    result.potential_V = std::move(potential);
    return result;
}

} // namespace vela
