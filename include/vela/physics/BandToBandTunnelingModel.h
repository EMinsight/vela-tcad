#pragma once

#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>
#include <string>

namespace vela {

/** Configuration for local band-to-band tunneling generation.
 *
 * The E2 defaults are the Sentaurus Device O-2018.06 silicon defaults.  The
 * stored units are always SI so the physical model is independent of Vela's
 * internal mesh/concentration unit system.
 */
struct BandToBandTunnelingConfig {
    std::string model = "none"; ///< "none" or the local Sentaurus-compatible "e2" model.
    Real prefactorA_SI = 3.4e23; ///< A [m^-1 s^-1 V^-2].
    Real exponentialB_V_per_m = 2.26e9; ///< B [V/m].
    Real minimumField_V_per_m = 0.0; ///< Hard evaluation floor; zero matches the E2 default.
    /// Spatial recovery/integration used to map the local E2 rate to nodes.
    /// The historical default evaluates one field per semiconductor cell.
    /// ``transport_node_lumped`` reconstructs a nodal field using only
    /// semiconductor-cell neighbours, then multiplies by the semiconductor
    /// lumped nodal area.
    std::string sourceIntegration = "semiconductor_cell_lumped";
    std::string jacobian = "frozen_field"; ///< Or exact but expensive "potential_finite_difference".
    Real jacobianRelativeStep = 1.0e-7;
};

class BandToBandTunnelingModel {
public:
    explicit BandToBandTunnelingModel(BandToBandTunnelingConfig config = {});

    bool enabled() const { return enabled_; }
    const BandToBandTunnelingConfig& config() const { return config_; }

    /// Local E2 pair-generation rate [m^-3 s^-1] for an SI electric field [V/m].
    Real generationRate(Real electricField_V_per_m) const;
    /// dG/d|F| in SI units, useful for independent derivative checks.
    Real generationRateDerivativeField(Real electricField_V_per_m) const;

private:
    BandToBandTunnelingConfig config_;
    bool enabled_ = false;
};

BandToBandTunnelingConfig bandToBandTunnelingConfigFromJson(
    const nlohmann::json& json,
    const char* context);

} // namespace vela
