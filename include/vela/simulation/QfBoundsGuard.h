#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/solver/GummelSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

enum class QfBoundsViolationMode {
    Warn,
    RejectAndRecover,
};

struct QfBoundsDiagnosticsConfig {
    bool enabled = true;
    QfBoundsViolationMode mode = QfBoundsViolationMode::Warn;
    Real margin_V = 0.5;
    // The quasi-Fermi potential is numerically unobservable when its
    // corresponding carrier density is effectively zero. Zero preserves the
    // historical all-node check; validation decks may select a physical floor.
    Real minCarrierDensity_m3 = 0.0;
    bool checkBandBending = false;
    Real builtInPotential_V = 0.0;
    std::string csvFile;
};

struct QfBoundsViolation {
    Index nodeId = 0;
    Real x = 0.0;
    Real y = 0.0;
    std::string variable;
    Real value = 0.0;
    Real lowerBound = 0.0;
    Real upperBound = 0.0;
    Real carrierDensity_m3 = std::numeric_limits<Real>::quiet_NaN();
};

struct QfBoundsEvaluation {
    bool checked = false;
    Real contactLower_V = 0.0;
    Real contactUpper_V = 0.0;
    Real margin_V = 0.0;
    Real bandBendingLimit_V = std::numeric_limits<Real>::infinity();
    std::vector<QfBoundsViolation> violations;

    bool valid() const { return violations.empty(); }
};

inline std::string toString(QfBoundsViolationMode mode)
{
    switch (mode) {
    case QfBoundsViolationMode::Warn:
        return "warn";
    case QfBoundsViolationMode::RejectAndRecover:
        return "reject_and_recover";
    }
    return "warn";
}

inline QfBoundsViolationMode qfBoundsViolationModeFromString(const std::string& mode)
{
    if (mode == "warn")
        return QfBoundsViolationMode::Warn;
    if (mode == "reject_and_recover")
        return QfBoundsViolationMode::RejectAndRecover;
    throw std::invalid_argument(
        "DCSweep: sweep.diagnostics.qf_bounds.mode must be 'warn' or 'reject_and_recover'.");
}

inline QfBoundsEvaluation evaluateQfBounds(
    const DeviceMesh& mesh,
    const DDSolution& solution,
    const std::unordered_map<std::string, Real>& contactBiases,
    const QfBoundsDiagnosticsConfig& config,
    Real activeBias_V)
{
    QfBoundsEvaluation eval;
    if (!config.enabled)
        return eval;

    eval.checked = true;
    eval.margin_V = config.margin_V;
    if (contactBiases.empty())
        return eval;

    Real minBias = std::numeric_limits<Real>::infinity();
    Real maxBias = -std::numeric_limits<Real>::infinity();
    for (const auto& [_, bias] : contactBiases) {
        if (!std::isfinite(bias))
            continue;
        minBias = std::min(minBias, bias);
        maxBias = std::max(maxBias, bias);
    }
    if (!std::isfinite(minBias) || !std::isfinite(maxBias))
        return eval;

    eval.contactLower_V = minBias - config.margin_V;
    eval.contactUpper_V = maxBias + config.margin_V;
    eval.bandBendingLimit_V =
        std::abs(activeBias_V) + std::abs(config.builtInPotential_V) + config.margin_V;

    const auto addViolation = [&](Index nodeId,
                                  std::string variable,
                                  Real value,
                                  Real lower,
                                  Real upper,
                                  Real carrierDensity_m3) {
        const Node& node = mesh.getNode(nodeId);
        eval.violations.push_back({
            nodeId,
            node.x,
            node.y,
            std::move(variable),
            value,
            lower,
            upper,
            carrierDensity_m3,
        });
    };

    const int nodeCount = static_cast<int>(mesh.numNodes());
    for (int i = 0; i < nodeCount; ++i) {
        const Index nodeId = static_cast<Index>(i);
        const Real phin = solution.phin(i);
        const Real phip = solution.phip(i);
        const Real electronDensity = i < solution.n.size()
            ? solution.n(i) : std::numeric_limits<Real>::quiet_NaN();
        const Real holeDensity = i < solution.p.size()
            ? solution.p(i) : std::numeric_limits<Real>::quiet_NaN();
        const bool checkElectronBounds =
            !std::isfinite(electronDensity) ||
            electronDensity > config.minCarrierDensity_m3;
        const bool checkHoleBounds =
            !std::isfinite(holeDensity) ||
            holeDensity > config.minCarrierDensity_m3;
        if (!std::isfinite(phin) ||
            (checkElectronBounds &&
             (phin < eval.contactLower_V || phin > eval.contactUpper_V))) {
            addViolation(
                nodeId, "phin", phin, eval.contactLower_V,
                eval.contactUpper_V, electronDensity);
        }
        if (!std::isfinite(phip) ||
            (checkHoleBounds &&
             (phip < eval.contactLower_V || phip > eval.contactUpper_V))) {
            addViolation(
                nodeId, "phip", phip, eval.contactLower_V,
                eval.contactUpper_V, holeDensity);
        }

        if (config.checkBandBending) {
            const Real psi = solution.psi(i);
            const Real psiMinusPhin = psi - phin;
            const Real phipMinusPsi = phip - psi;
            const Real lower = -eval.bandBendingLimit_V;
            const Real upper = eval.bandBendingLimit_V;
            if (!std::isfinite(psiMinusPhin) ||
                psiMinusPhin < lower ||
                psiMinusPhin > upper) {
                addViolation(
                    nodeId, "psi_minus_phin", psiMinusPhin, lower, upper,
                    electronDensity);
            }
            if (!std::isfinite(phipMinusPsi) ||
                phipMinusPsi < lower ||
                phipMinusPsi > upper) {
                addViolation(
                    nodeId, "phip_minus_psi", phipMinusPsi, lower, upper,
                    holeDensity);
            }
        }
    }

    return eval;
}

inline DDSolution resetQfBoundsViolationsToNearestContactBias(
    const DDSolution& solution,
    const QfBoundsEvaluation& eval)
{
    DDSolution reset = solution;
    for (const QfBoundsViolation& violation : eval.violations) {
        const Real target =
            std::abs(violation.value - eval.contactLower_V) <
                std::abs(violation.value - eval.contactUpper_V)
            ? eval.contactLower_V + eval.margin_V
            : eval.contactUpper_V - eval.margin_V;
        const int i = static_cast<int>(violation.nodeId);
        if (violation.variable == "phin" && i < reset.phin.size())
            reset.phin(i) = target;
        else if (violation.variable == "phip" && i < reset.phip.size())
            reset.phip(i) = target;
    }
    return reset;
}

} // namespace vela
