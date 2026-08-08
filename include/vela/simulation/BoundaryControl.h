#pragma once

#include "vela/core/Types.h"

#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace vela::detail {

struct MonotoneBoundaryRootBracket {
    Real negativeVoltage = 0.0;
    Real negativeResidual = 0.0;
    Real positiveVoltage = 0.0;
    Real positiveResidual = 0.0;
};

struct MonotoneBoundaryRootConfig {
    Real maxStep = 0.1;
    Real predictorMaxStepFactor = 4.0;
    Real residualTolerance = 1.0e-8;
    Real voltageTolerance = 1.0e-8;
    int maxBracketSteps = 200;
    int maxIterations = 40;
    /// Residual at initialVoltage when the caller already owns the matching
    /// converged device state. Reusing it avoids an identical expensive solve.
    std::optional<Real> initialResidual;
    /// Preferred first trial, normally predicted from the previous two
    /// accepted boundary-control operating points.
    std::optional<Real> predictedVoltage;
    /// Previously evaluated sign-changing bracket. Supplying it resumes the
    /// scalar corrector without rerunning either expensive device solve.
    std::optional<MonotoneBoundaryRootBracket> initialBracket;
};

struct MonotoneBoundaryRootResult {
    Real voltage = 0.0;
    Real residual = 0.0;
    int evaluations = 0;
    bool converged = false;
};

using BoundaryResidualEvaluator = std::function<Real(Real)>;

struct BoundaryVoltagePrediction {
    Real voltage = 0.0;
    bool curvatureAccelerated = false;
};

/// Predict the next boundary voltage from accepted (target, voltage) pairs.
/// Two points give a secant predictor. Four equally spaced targets additionally
/// support a guarded second-order trend in the shrinking voltage increments.
std::optional<BoundaryVoltagePrediction> predictBoundaryVoltageFromHistory(
    const std::vector<std::pair<Real, Real>>& history,
    Real target);

/// Solve a scalar boundary equation whose residual increases with device
/// voltage. The evaluator may run a complete device solve at each voltage.
MonotoneBoundaryRootResult solveMonotoneBoundaryRoot(
    Real initialVoltage,
    const MonotoneBoundaryRootConfig& config,
    const BoundaryResidualEvaluator& evaluate);

/// Sentaurus-compatible 2-D series-resistor convention. Resistance is in
/// ohm*um and current is in A/um, so the product is a voltage drop.
Real externalResistorOuterVoltage(Real innerVoltage_V,
                                  Real resistance_ohm_um,
                                  Real directedCurrent_A_per_um);

Real externalResistorLoadLineResidual(Real innerVoltage_V,
                                      Real outerVoltage_V,
                                      Real resistance_ohm_um,
                                      Real directedCurrent_A_per_um);

} // namespace vela::detail
