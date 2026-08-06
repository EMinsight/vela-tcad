#include "vela/physics/BandToBandTunnelingModel.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace vela {

BandToBandTunnelingModel::BandToBandTunnelingModel(
    BandToBandTunnelingConfig config)
    : config_(std::move(config))
{
    if (config_.model == "none") {
        enabled_ = false;
    } else if (config_.model == "e2" || config_.model == "sentaurus_e2") {
        enabled_ = true;
        config_.model = "e2";
    } else {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: model must be 'none' or 'e2'.");
    }

    if (!(config_.prefactorA_SI >= 0.0) || !std::isfinite(config_.prefactorA_SI))
        throw std::invalid_argument(
            "BandToBandTunnelingModel: prefactor A must be finite and non-negative.");
    if (!(config_.exponentialB_V_per_m > 0.0) ||
        !std::isfinite(config_.exponentialB_V_per_m)) {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: exponential B must be positive and finite.");
    }
    if (!(config_.minimumField_V_per_m >= 0.0) ||
        !std::isfinite(config_.minimumField_V_per_m)) {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: minimum field must be finite and non-negative.");
    }
    if (config_.sourceIntegration != "semiconductor_cell_lumped" &&
        config_.sourceIntegration != "transport_node_lumped") {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: source_integration must be "
            "'semiconductor_cell_lumped' or 'transport_node_lumped'.");
    }
    if (config_.jacobian != "potential_finite_difference" &&
        config_.jacobian != "frozen_field") {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: jacobian must be "
            "'potential_finite_difference' or 'frozen_field'.");
    }
    if (!(config_.jacobianRelativeStep > 0.0) ||
        !std::isfinite(config_.jacobianRelativeStep)) {
        throw std::invalid_argument(
            "BandToBandTunnelingModel: Jacobian relative step must be positive and finite.");
    }
}

Real BandToBandTunnelingModel::generationRate(Real electricField_V_per_m) const
{
    if (!enabled_ || !std::isfinite(electricField_V_per_m))
        return 0.0;
    const Real field = std::abs(electricField_V_per_m);
    if (field <= config_.minimumField_V_per_m || field == 0.0)
        return 0.0;

    const Real exponent = -config_.exponentialB_V_per_m / field;
    if (exponent < -745.0)
        return 0.0;
    return config_.prefactorA_SI * field * field * std::exp(exponent);
}

Real BandToBandTunnelingModel::generationRateDerivativeField(
    Real electricField_V_per_m) const
{
    if (!enabled_ || !std::isfinite(electricField_V_per_m))
        return 0.0;
    const Real field = std::abs(electricField_V_per_m);
    if (field <= config_.minimumField_V_per_m || field == 0.0)
        return 0.0;
    const Real rate = generationRate(field);
    return rate * (2.0 / field + config_.exponentialB_V_per_m / (field * field));
}

BandToBandTunnelingConfig bandToBandTunnelingConfigFromJson(
    const nlohmann::json& json,
    const char* context)
{
    BandToBandTunnelingConfig config;
    if (json.is_string()) {
        config.model = json.get<std::string>();
        return config;
    }
    if (!json.is_object()) {
        throw std::invalid_argument(
            std::string(context) + ": band_to_band must be a string or object.");
    }

    config.model = json.value("model", config.model);
    config.prefactorA_SI = json.value(
        "A_m_inv_s_inv_V_inv2", config.prefactorA_SI);
    // Sentaurus-native aliases make direct transcription of parameter files
    // explicit while retaining SI keys as the canonical Vela interface.
    if (json.contains("A_cm_inv_s_inv_V_inv2")) {
        config.prefactorA_SI =
            100.0 * json.at("A_cm_inv_s_inv_V_inv2").get<Real>();
    }
    config.exponentialB_V_per_m = json.value(
        "B_V_per_m", config.exponentialB_V_per_m);
    if (json.contains("B_V_per_cm")) {
        config.exponentialB_V_per_m =
            100.0 * json.at("B_V_per_cm").get<Real>();
    }
    config.minimumField_V_per_m = json.value(
        "minimum_field_V_per_m", config.minimumField_V_per_m);
    if (json.contains("minimum_field_V_per_cm")) {
        config.minimumField_V_per_m =
            100.0 * json.at("minimum_field_V_per_cm").get<Real>();
    }
    config.sourceIntegration = json.value(
        "source_integration", config.sourceIntegration);
    config.jacobian = json.value("jacobian", config.jacobian);
    config.jacobianRelativeStep = json.value(
        "jacobian_relative_step", config.jacobianRelativeStep);
    return config;
}

} // namespace vela
