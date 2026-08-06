#include "vela/post/ContactCurrent.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vela {
namespace {

Real validatedThermalVoltage(Real temperature_K)
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("ContactCurrent: temperature_K must be positive.");
    return constants::kb * temperature_K / constants::q;
}

} // namespace

ContactCurrent::ContactCurrent(const DeviceMesh& mesh,
                               const MaterialDatabase& matdb,
                               const DopingModel& doping,
                               MobilityModelConfig mobilityConfig,
                               Real temperature_K,
                               DDScalingSpec scaling,
                               BandgapNarrowingConfig bandgapNarrowingConfig,
                               CarrierStatisticsConfig carrierStatistics)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , couple_(detail::computeEdgeCouplings(mesh))
    , mobilityConfig_(mobilityConfig)
    , mobility_(makeMobilityModel(mobilityConfig))
    , thermalVoltage_(validatedThermalVoltage(temperature_K))
    , scaling_(scaling)
    , ni_(detail::buildValidatedEffectiveNodeNi(
          "ContactCurrent",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          validatedThermalVoltage(temperature_K)))
    , Nc_(detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, true))
    , Nv_(detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, false))
    , carrierStatistics_(std::move(carrierStatistics))
{}


ContactCurrentResult ContactCurrent::compute(const DDSolution& solution,
                                             const std::string& contactName) const
{
    return computeDetailed(solution, contactName).totals;
}

ContactCurrentResult ContactCurrent::compute(
    const DDSolution& solution,
    const std::string& contactName,
    const ContactCurrentEdgeOverrides& overrides) const
{
    return computeDetailed(solution, contactName, overrides).totals;
}

ContactCurrentDetailedResult ContactCurrent::computeDetailed(
    const DDSolution& solution,
    const std::string& contactName) const
{
    static const ContactCurrentEdgeOverrides noOverrides;
    return computeDetailed(solution, contactName, noOverrides);
}

ContactCurrentResult ContactCurrent::computeFromResidual(
    const CoupledDDAssembler& assembler,
    const VectorXd& x,
    const std::string& contactName) const
{
    const Contact* contact = nullptr;
    for (const Contact& candidate : mesh_.contacts()) {
        if (candidate.name == contactName) {
            contact = &candidate;
            break;
        }
    }
    if (contact == nullptr)
        throw std::invalid_argument("ContactCurrent: unknown contact '" + contactName + "'.");

    const int N = static_cast<int>(assembler.numNodes());
    const int phinOffset = N;
    const int phipOffset = 2 * N;
    const VectorXd r = assembler.residual(x, CoupledDDBoundaryConditions{});
    const Real continuityScale = assembler.continuityResidualScale();
    const Real currentLineFactor = scaling_.enabled
        ? scaling_.currentDensityLineIntegralFactor
        : 1.0;

    ContactCurrentResult result;
    for (Index node : contact->node_ids) {
        const int row = static_cast<int>(node);
        result.electronCurrent -= constants::q * r(phinOffset + row) * continuityScale * currentLineFactor;
        result.holeCurrent -= constants::q * r(phipOffset + row) * continuityScale * currentLineFactor;
    }
    result.totalCurrent = result.electronCurrent - result.holeCurrent;
    return result;
}
ContactCurrentDetailedResult ContactCurrent::computeDetailed(
    const DDSolution& solution,
    const std::string& contactName,
    const ContactCurrentEdgeOverrides& overrides) const
{
    const Contact* contact = nullptr;
    for (const Contact& candidate : mesh_.contacts()) {
        if (candidate.name == contactName) {
            contact = &candidate;
            break;
        }
    }
    if (contact == nullptr)
        throw std::invalid_argument("ContactCurrent: unknown contact '" + contactName + "'.");

    std::unordered_set<Index> contactNodes(contact->node_ids.begin(), contact->node_ids.end());
    const Real temperature_K = thermalVoltage_ * constants::q / constants::kb;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh_, matdb_, temperature_K);
    const Real fieldFactor = scaling_.enabled
        ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const bool hasReferencedElectronQf =
        solution.phinIncrement.size() == solution.phin.size();
    const bool hasReferencedHoleQf =
        solution.phipIncrement.size() == solution.phip.size();
    const VectorXd& electronQf = hasReferencedElectronQf
        ? solution.phinIncrement : solution.phin;
    const VectorXd& holeQf = hasReferencedHoleQf
        ? solution.phipIncrement : solution.phip;
    const bool vectorQfMobility =
        mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient" &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, electronQf, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials, holeQf, fieldFactor)
        : std::vector<Real>{};

    ContactCurrentDetailedResult detailed;
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const bool n0OnContact = contactNodes.count(edge.n0) > 0;
        const bool n1OnContact = contactNodes.count(edge.n1) > 0;
        if (n0OnContact == n1OnContact)
            continue;
        if (edge.length < 1.0e-30 || couple_[e] <= 0.0)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);

        // The solver returns physical DDSolution fields in the active internal unit system.
        const Real psi_i = solution.psi(i);
        const Real psi_j = solution.psi(j);
        const Real n_i = solution.n(i);
        const Real n_j = solution.n(j);
        const Real p_i = solution.p(i);
        const Real p_j = solution.p(j);
        const Real dpsi = psi_j - psi_i;
        const Real edgeLength = edge.length;

        const Real currentLineFactor = scaling_.enabled
            ? scaling_.currentDensityLineIntegralFactor
            : 1.0;
        const Real electricField = std::abs(dpsi / edgeLength) * fieldFactor;
        const Real phin_i = hasReferencedElectronQf
            ? solution.phinIncrement(i) : solution.phin(i);
        const Real phin_j = hasReferencedElectronQf
            ? solution.phinIncrement(j) : solution.phin(j);
        const Real phip_i = hasReferencedHoleQf
            ? solution.phipIncrement(i) : solution.phip(i);
        const Real phip_j = hasReferencedHoleQf
            ? solution.phipIncrement(j) : solution.phip(j);
        const Real electronPsi_i = psi_i -
            (hasReferencedElectronQf ? solution.electronQfReference_V : 0.0);
        const Real electronPsi_j = psi_j -
            (hasReferencedElectronQf ? solution.electronQfReference_V : 0.0);
        const Real holePsi_i = psi_i -
            (hasReferencedHoleQf ? solution.holeQfReference_V : 0.0);
        const Real holePsi_j = psi_j -
            (hasReferencedHoleQf ? solution.holeQfReference_V : 0.0);
        Real phip_i_forHole = phip_i;
        Real phip_j_forHole = phip_j;
        bool holeQfDropOverrideApplied = false;
        const auto holeDropIt = overrides.holeQuasiFermiDropByEdge.find(e);
        if (holeDropIt != overrides.holeQuasiFermiDropByEdge.end()
            && std::isfinite(holeDropIt->second)) {
            phip_j_forHole = phip_i_forHole + holeDropIt->second;
            holeQfDropOverrideApplied = true;
        }
        const Real electronMobilityField =
            vectorQfMobility
            ? electronVectorMobilityFields[e]
            : mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient"
            ? std::abs((phin_j - phin_i) / edgeLength) * fieldFactor
            : electricField;
        const Real holeMobilityField =
            vectorQfMobility
            ? holeVectorMobilityFields[e]
            : mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient"
            ? std::abs((phip_j - phip_i) / edgeLength) * fieldFactor
            : electricField;

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &solution.psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &solution.psi);

        // SG fluxes in physical units.  Mirror CoupledDDAssembler residual:
        // use the cancellation-free quasi-Fermi balanced form, including the
        // variable-ni generalization needed for BGN/effective-ni edges.  The
        // density-based form B(-u)*n0 - B(+u)*n1 does not cancel flat
        // quasi-Fermi levels when ni varies across the edge.
        const Index idxI = edge.n0;
        const Index idxJ = edge.n1;
        const Real ni_i = ni_[idxI];
        const Real ni_j = ni_[idxJ];
        const bool fermiDirac = usesFermiDirac(carrierStatistics_);
        Real electronGeneralizedFactor = 1.0;
        Real electronGeneralizedArgument = dpsi / thermalVoltage_;
        Real holeGeneralizedFactor = 1.0;
        Real holeGeneralizedArgument = dpsi / thermalVoltage_;
        Real electronContinuityFlux01 = 0.0;
        Real electronFlux01 = 0.0;
        if (mun > 0.0) {
            const Real coef = mun * thermalVoltage_ * fieldFactor / edgeLength;
            if (fermiDirac) {
                const Real etaI = (electronPsi_i - phin_i) / thermalVoltage_
                    + std::log(ni_i / Nc_[idxI]);
                const Real etaJ = (electronPsi_j - phin_j) / thermalVoltage_
                    + std::log(ni_j / Nc_[idxJ]);
                const Real driftPotential = electronPsi_j - electronPsi_i
                    + thermalVoltage_ * std::log(
                        (ni_j / Nc_[idxJ]) / (ni_i / Nc_[idxI]));
                electronGeneralizedFactor = sgGeneralizedEinsteinFactor(
                    n_i, n_j, etaI, etaJ);
                electronGeneralizedArgument = driftPotential /
                    (thermalVoltage_ * electronGeneralizedFactor);
                electronContinuityFlux01 = sgElectronFermiDiracContinuityFlux(
                    n_i, n_j, etaI, etaJ, driftPotential,
                    phin_i, phin_j, thermalVoltage_, coef);
            } else {
                electronContinuityFlux01 = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_i, ni_j, electronPsi_i, electronPsi_j,
                    phin_i, phin_j, thermalVoltage_, coef);
            }
            // sgElectronFlux = -sgElectronContinuityFlux by definition.
            electronFlux01 = -electronContinuityFlux01;
        }
        Real holeContinuityFlux01 = 0.0;
        Real holeFlux01 = 0.0;
        if (mup > 0.0) {
            const Real coef = mup * thermalVoltage_ * fieldFactor / edgeLength;
            if (fermiDirac) {
                const Real etaI = (phip_i_forHole - holePsi_i) / thermalVoltage_
                    + std::log(ni_i / Nv_[idxI]);
                const Real etaJ = (phip_j_forHole - holePsi_j) / thermalVoltage_
                    + std::log(ni_j / Nv_[idxJ]);
                const Real driftPotential = holePsi_j - holePsi_i
                    + thermalVoltage_ * std::log(
                        (ni_i / Nv_[idxI]) / (ni_j / Nv_[idxJ]));
                holeGeneralizedFactor = sgGeneralizedEinsteinFactor(
                    p_i, p_j, etaI, etaJ);
                holeGeneralizedArgument = driftPotential /
                    (thermalVoltage_ * holeGeneralizedFactor);
                holeContinuityFlux01 = sgHoleFermiDiracContinuityFlux(
                    p_i, p_j, etaI, etaJ, driftPotential,
                    phip_i_forHole, phip_j_forHole,
                    thermalVoltage_, coef);
            } else {
                holeContinuityFlux01 = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_i, ni_j, holePsi_i, holePsi_j,
                    phip_i_forHole, phip_j_forHole,
                    thermalVoltage_, coef);
            }
            holeFlux01 = -holeContinuityFlux01;
        }

        // Algebraic SG split: J = J_drift + J_diffusion.
        const SGEdgeWeights weights = sgEdgeWeights(dpsi, thermalVoltage_);
        const Real bAvg = 0.5 * (weights.b_plus + weights.b_minus);
        Real electronDriftFlux01 = 0.0;
        Real electronDiffusionFlux01 = 0.0;
        Real holeDriftFlux01 = 0.0;
        Real holeDiffusionFlux01 = 0.0;
        if (fermiDirac) {
            if (mun > 0.0) {
                const Real bMinus = bernoulli(-electronGeneralizedArgument);
                const Real bPlus = bernoulli(electronGeneralizedArgument);
                const Real coefficient = mun * thermalVoltage_ * fieldFactor /
                    edgeLength * electronGeneralizedFactor;
                electronDriftFlux01 = -coefficient * 0.5 * (bMinus - bPlus)
                    * (n_i + n_j);
                electronDiffusionFlux01 = -coefficient * 0.5 * (bMinus + bPlus)
                    * (n_i - n_j);
            }
            if (mup > 0.0) {
                const Real bPlus = bernoulli(holeGeneralizedArgument);
                const Real bMinus = bernoulli(-holeGeneralizedArgument);
                const Real coefficient = mup * thermalVoltage_ * fieldFactor /
                    edgeLength * holeGeneralizedFactor;
                holeDriftFlux01 = -coefficient * 0.5 * (bPlus - bMinus)
                    * (p_i + p_j);
                holeDiffusionFlux01 = -coefficient * 0.5 * (bPlus + bMinus)
                    * (p_i - p_j);
            }
        } else {
            electronDriftFlux01 = (mun > 0.0)
                ? mun * (dpsi / edgeLength) * fieldFactor * (0.5 * (n_i + n_j))
                : 0.0;
            electronDiffusionFlux01 = (mun > 0.0)
                ? mun * (thermalVoltage_ / edgeLength) * fieldFactor * bAvg * (n_i - n_j)
                : 0.0;
            holeDriftFlux01 = (mup > 0.0)
                ? mup * (dpsi / edgeLength) * fieldFactor * (0.5 * (p_i + p_j))
                : 0.0;
            holeDiffusionFlux01 = (mup > 0.0)
                ? mup * (thermalVoltage_ / edgeLength) * fieldFactor * bAvg * (p_j - p_i)
                : 0.0;
        }

        const Real outwardSign = n0OnContact ? 1.0 : -1.0;
        // Current density in the active internal unit system times edge length gives current per internal device depth.
        const Real electronCurrent = constants::q * outwardSign * electronFlux01 * couple_[e] * currentLineFactor;
        const Real electronDriftCurrent = constants::q * outwardSign * electronDriftFlux01 * couple_[e] * currentLineFactor;
        const Real electronDiffusionCurrent = constants::q * outwardSign * electronDiffusionFlux01 * couple_[e] * currentLineFactor;
        const Real holeCurrent = constants::q * outwardSign * holeFlux01 * couple_[e] * currentLineFactor;
        const Real holeDriftCurrent = constants::q * outwardSign * holeDriftFlux01 * couple_[e] * currentLineFactor;
        const Real holeDiffusionCurrent = constants::q * outwardSign * holeDiffusionFlux01 * couple_[e] * currentLineFactor;

        detailed.totals.electronCurrent += electronCurrent;
        detailed.totals.electronDriftCurrent += electronDriftCurrent;
        detailed.totals.electronDiffusionCurrent += electronDiffusionCurrent;
        detailed.totals.holeCurrent += holeCurrent;
        detailed.totals.holeDriftCurrent += holeDriftCurrent;
        detailed.totals.holeDiffusionCurrent += holeDiffusionCurrent;

        ContactCurrentEdgeDiagnostic edgeDiag;
        edgeDiag.edgeId = e;
        edgeDiag.node0 = edge.n0;
        edgeDiag.node1 = edge.n1;
        edgeDiag.edgeLength_m = scaling_.unitSystem.internalLengthToMeters(edge.length);
        edgeDiag.edgeCouple_m = scaling_.unitSystem.internalLengthToMeters(couple_[e]);
        edgeDiag.outwardSign = outwardSign;
        edgeDiag.bernoulliU = dpsi / thermalVoltage_;
        edgeDiag.bernoulliBplus = weights.b_plus;
        edgeDiag.bernoulliBminus = weights.b_minus;
        edgeDiag.electronUsedQuasiFermi = true;
        edgeDiag.holeUsedQuasiFermi = true;
        edgeDiag.psi0 = psi_i;
        edgeDiag.psi1 = psi_j;
        edgeDiag.phin0 = phin_i +
            (hasReferencedElectronQf ? solution.electronQfReference_V : 0.0);
        edgeDiag.phin1 = phin_j +
            (hasReferencedElectronQf ? solution.electronQfReference_V : 0.0);
        edgeDiag.phip0 = phip_i_forHole +
            (hasReferencedHoleQf ? solution.holeQfReference_V : 0.0);
        edgeDiag.phip1 = phip_j_forHole +
            (hasReferencedHoleQf ? solution.holeQfReference_V : 0.0);
        edgeDiag.holeQfDropOverrideApplied = holeQfDropOverrideApplied;
        edgeDiag.n0 = n_i;
        edgeDiag.n1 = n_j;
        edgeDiag.p0 = p_i;
        edgeDiag.p1 = p_j;
        edgeDiag.ni0 = ni_i;
        edgeDiag.ni1 = ni_j;
        edgeDiag.mun = mun;
        edgeDiag.mup = mup;
        edgeDiag.electronContinuityFlux = electronContinuityFlux01;
        edgeDiag.holeContinuityFlux = holeContinuityFlux01;
        edgeDiag.electronCurrent = electronCurrent;
        edgeDiag.electronDriftCurrent = electronDriftCurrent;
        edgeDiag.electronDiffusionCurrent = electronDiffusionCurrent;
        edgeDiag.holeCurrent = holeCurrent;
        edgeDiag.holeDriftCurrent = holeDriftCurrent;
        edgeDiag.holeDiffusionCurrent = holeDiffusionCurrent;
        edgeDiag.totalCurrent = electronCurrent - holeCurrent;
        detailed.edges.push_back(std::move(edgeDiag));
    }

    // Sign convention: electronCurrent and holeCurrent accumulate
    //   q * (carrier-particle inflow into the contact from the device) * couple.
    // With the electron carrier charge being -q, the contribution of electrons
    // to the conventional current supplied into the contact from the external
    // circuit is -(particle inflow), so the total terminal current is
    //   I_total = I_electron - I_hole.
    // Using `I_electron + I_hole` (as previously) double-adds the volume
    // recombination integral into the terminal current and breaks the
    // Kirchhoff balance |I_anode| = |I_cathode| for a two-terminal device.
    detailed.totals.totalCurrent = detailed.totals.electronCurrent - detailed.totals.holeCurrent;
    return detailed;
}


ContactCurrentResult ContactCurrent::compute(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const DDSolution& solution,
    const std::string& contactName,
    const MobilityModelConfig& mobilityConfig,
    Real temperature_K,
    DDScalingSpec scaling,
    const BandgapNarrowingConfig& bandgapNarrowingConfig,
    const CarrierStatisticsConfig& carrierStatistics)
{
    return ContactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K, scaling,
                          bandgapNarrowingConfig, carrierStatistics)
        .compute(solution, contactName);
}

} // namespace vela
