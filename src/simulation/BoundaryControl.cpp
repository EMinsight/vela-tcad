#include "vela/simulation/BoundaryControl.h"
#include "vela/core/PerformanceProfiler.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela::detail {

namespace {

void validate(const MonotoneBoundaryRootConfig& config)
{
    if (!std::isfinite(config.maxStep) || config.maxStep <= 0.0)
        throw std::invalid_argument("Boundary control: maxStep must be finite and positive.");
    if (!std::isfinite(config.predictorMaxStepFactor) ||
        config.predictorMaxStepFactor < 1.0) {
        throw std::invalid_argument(
            "Boundary control: predictorMaxStepFactor must be finite and at least one.");
    }
    if (!std::isfinite(config.residualTolerance) || config.residualTolerance <= 0.0)
        throw std::invalid_argument(
            "Boundary control: residualTolerance must be finite and positive.");
    if (!std::isfinite(config.voltageTolerance) || config.voltageTolerance <= 0.0)
        throw std::invalid_argument(
            "Boundary control: voltageTolerance must be finite and positive.");
    if (config.maxBracketSteps <= 0)
        throw std::invalid_argument("Boundary control: maxBracketSteps must be positive.");
    if (config.maxIterations <= 0)
        throw std::invalid_argument("Boundary control: maxIterations must be positive.");
}

} // namespace

std::optional<BoundaryVoltagePrediction> predictBoundaryVoltageFromHistory(
    const std::vector<std::pair<Real, Real>>& history,
    Real target)
{
    if (history.size() < 2)
        return std::nullopt;
    const auto& [target0, voltage0] = history.at(history.size() - 2);
    const auto& [target1, voltage1] = history.back();
    const Real targetStep = target1 - target0;
    if (!std::isfinite(targetStep) || targetStep == 0.0)
        return std::nullopt;

    BoundaryVoltagePrediction prediction;
    prediction.voltage = voltage1 +
        (target - target1) * (voltage1 - voltage0) / targetStep;
    if (history.size() < 4)
        return prediction;

    const auto& [targetM3, voltageM3] = history.at(history.size() - 4);
    const auto& [targetM2, voltageM2] = history.at(history.size() - 3);
    const Real stepM2 = targetM2 - targetM3;
    const Real stepM1 = target0 - targetM2;
    const Real scale = std::max({
        std::abs(stepM2), std::abs(stepM1), std::abs(targetStep),
        std::abs(target - target1), Real{1.0}});
    const Real tolerance = 1.0e-10 * scale;
    if (std::abs(stepM2 - targetStep) > tolerance ||
        std::abs(stepM1 - targetStep) > tolerance ||
        std::abs((target - target1) - targetStep) > tolerance) {
        return prediction;
    }

    const Real delta0 = voltageM2 - voltageM3;
    const Real delta1 = voltage0 - voltageM2;
    const Real delta2 = voltage1 - voltage0;
    if (!std::isfinite(delta0) || !std::isfinite(delta1) ||
        !std::isfinite(delta2) || delta0 == 0.0 || delta1 == 0.0 ||
        std::signbit(delta0) != std::signbit(delta1) ||
        std::signbit(delta1) != std::signbit(delta2)) {
        return prediction;
    }
    const Real ratio0 = delta1 / delta0;
    const Real ratio1 = delta2 / delta1;
    const Real predictedRatio = 2.0 * ratio1 - ratio0;
    if (!std::isfinite(predictedRatio) ||
        predictedRatio < 0.25 || predictedRatio > 1.25) {
        return prediction;
    }
    const Real accelerated = voltage1 + delta2 * predictedRatio;
    const Real linearDelta = prediction.voltage - voltage1;
    const Real acceleratedDelta = accelerated - voltage1;
    if (!std::isfinite(accelerated) ||
        std::signbit(acceleratedDelta) != std::signbit(linearDelta) ||
        std::abs(acceleratedDelta) < 0.25 * std::abs(linearDelta) ||
        std::abs(acceleratedDelta) > 1.25 * std::abs(linearDelta)) {
        return prediction;
    }
    prediction.voltage = accelerated;
    prediction.curvatureAccelerated = true;
    return prediction;
}

MonotoneBoundaryRootResult solveMonotoneBoundaryRoot(
    Real initialVoltage,
    const MonotoneBoundaryRootConfig& config,
    const BoundaryResidualEvaluator& evaluate)
{
    ScopedPerformanceTimer timer("boundary.root_total");
    incrementPerformanceCounter("boundary.root_calls");
    validate(config);
    if (!std::isfinite(initialVoltage))
        throw std::invalid_argument("Boundary control: initial voltage must be finite.");

    MonotoneBoundaryRootResult result;
    auto evaluated = [&](Real voltage) {
        incrementPerformanceCounter("boundary.residual_evaluations");
        const Real residual = evaluate(voltage);
        ++result.evaluations;
        if (!std::isfinite(residual))
            throw std::runtime_error("Boundary control: residual evaluator returned a non-finite value.");
        return residual;
    };

    Real lowerX = 0.0;
    Real lowerF = 0.0;
    Real upperX = 0.0;
    Real upperF = 0.0;
    bool bracketed = false;
    if (config.initialBracket.has_value()) {
        const auto& bracket = *config.initialBracket;
        if (!std::isfinite(bracket.negativeVoltage) ||
            !std::isfinite(bracket.negativeResidual) ||
            !std::isfinite(bracket.positiveVoltage) ||
            !std::isfinite(bracket.positiveResidual) ||
            bracket.negativeResidual > 0.0 ||
            bracket.positiveResidual < 0.0 ||
            bracket.negativeVoltage >= bracket.positiveVoltage) {
            throw std::invalid_argument(
                "Boundary control: initialBracket must contain finite, ordered, "
                "sign-changing endpoints.");
        }
        lowerX = bracket.negativeVoltage;
        lowerF = bracket.negativeResidual;
        upperX = bracket.positiveVoltage;
        upperF = bracket.positiveResidual;
        if (std::abs(lowerF) <= config.residualTolerance) {
            result.voltage = lowerX;
            result.residual = lowerF;
            result.converged = true;
            return result;
        }
        if (std::abs(upperF) <= config.residualTolerance) {
            result.voltage = upperX;
            result.residual = upperF;
            result.converged = true;
            return result;
        }
        bracketed = true;
    }

    Real x0 = initialVoltage;
    Real f0 = 0.0;
    if (!bracketed) {
        f0 = config.initialResidual.has_value()
            ? *config.initialResidual
            : evaluated(x0);
        if (!std::isfinite(f0))
            throw std::invalid_argument(
                "Boundary control: initialResidual must be finite when specified.");
    }
    if (!bracketed && std::abs(f0) <= config.residualTolerance) {
        result.voltage = x0;
        result.residual = f0;
        result.converged = true;
        return result;
    }

    Real x1 = x0;
    Real f1 = f0;
    if (!bracketed) {
        Real direction = f0 < 0.0 ? 1.0 : -1.0;
        x1 = x0 + direction * config.maxStep;
        if (config.predictedVoltage.has_value() &&
            std::isfinite(*config.predictedVoltage) &&
            direction * (*config.predictedVoltage - x0) > 0.0) {
            const Real predictedDelta = std::clamp(
                *config.predictedVoltage - x0,
                -config.maxStep * config.predictorMaxStepFactor,
                config.maxStep * config.predictorMaxStepFactor);
            x1 = x0 + predictedDelta;
        }
        f1 = evaluated(x1);
        if (std::abs(f1) <= config.residualTolerance) {
            result.voltage = x1;
            result.residual = f1;
            result.converged = true;
            return result;
        }
        for (int step = 0; step < config.maxBracketSteps; ++step) {
            incrementPerformanceCounter("boundary.bracket_steps");
            if ((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)) {
                bracketed = true;
                break;
            }
            direction = f1 < 0.0 ? 1.0 : -1.0;
            Real next = std::numeric_limits<Real>::quiet_NaN();
            const Real denominator = f1 - f0;
            if (std::isfinite(denominator) && denominator != 0.0)
                next = x1 - f1 * (x1 - x0) / denominator;
            const Real maxPredictorStep =
                config.maxStep * config.predictorMaxStepFactor;
            if (!std::isfinite(next) || direction * (next - x1) <= 0.0) {
                next = x1 + direction * config.maxStep;
            } else {
                next = x1 + std::clamp(
                    next - x1, -maxPredictorStep, maxPredictorStep);
            }
            x0 = x1;
            f0 = f1;
            x1 = next;
            f1 = evaluated(x1);
            if (std::abs(f1) <= config.residualTolerance) {
                result.voltage = x1;
                result.residual = f1;
                result.converged = true;
                return result;
            }
        }
    }
    if (!bracketed)
        throw std::runtime_error("Boundary control: failed to bracket the scalar boundary root.");

    if (!config.initialBracket.has_value()) {
        lowerX = x0;
        lowerF = f0;
        upperX = x1;
        upperF = f1;
        if (lowerX > upperX) {
            std::swap(lowerX, upperX);
            std::swap(lowerF, upperF);
        }
        if (lowerF > 0.0) {
            std::swap(lowerX, upperX);
            std::swap(lowerF, upperF);
        }
    }

    for (int iteration = 0; iteration < config.maxIterations; ++iteration) {
        incrementPerformanceCounter("boundary.root_updates");
        const Real width = upperX - lowerX;
        Real trial = lowerX - lowerF * width / (upperF - lowerF);
        const Real guard = std::max(1.0e-12, 1.0e-6 * width);
        if (!std::isfinite(trial)) {
            trial = 0.5 * (lowerX + upperX);
        } else {
            trial = std::clamp(trial, lowerX + guard, upperX - guard);
        }
        const Real residual = evaluated(trial);
        result.voltage = trial;
        result.residual = residual;
        if (std::abs(residual) <= config.residualTolerance ||
            width <= config.voltageTolerance) {
            result.converged = true;
            return result;
        }
        if (residual < 0.0) {
            lowerX = trial;
            lowerF = residual;
        } else {
            upperX = trial;
            upperF = residual;
        }
    }

    throw std::runtime_error("Boundary control: scalar boundary root exceeded iteration limit.");
}

Real externalResistorOuterVoltage(Real innerVoltage_V,
                                  Real resistance_ohm_um,
                                  Real directedCurrent_A_per_um)
{
    return innerVoltage_V + resistance_ohm_um * directedCurrent_A_per_um;
}

Real externalResistorLoadLineResidual(Real innerVoltage_V,
                                      Real outerVoltage_V,
                                      Real resistance_ohm_um,
                                      Real directedCurrent_A_per_um)
{
    return externalResistorOuterVoltage(
               innerVoltage_V, resistance_ohm_um, directedCurrent_A_per_um) -
        outerVoltage_V;
}

} // namespace vela::detail
