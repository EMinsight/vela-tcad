#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/BVProcessProbe.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/Bernoulli.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/equation/AssemblerUtils.h"
#include "vela/physics/CarrierStatistics.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace vela {
namespace {

Real bernoulliDerivative(Real x)
{
    const Real ax = std::abs(x);
    if (ax < 1.0e-8)
        return -0.5 + x / 6.0 - x * x * x / 180.0;
    if (x > 500.0)
        return (1.0 - x) * std::exp(-x);
    if (x < -500.0)
        return -1.0;

    const Real em1 = std::expm1(x);
    const Real ex = em1 + 1.0;
    return (em1 - x * ex) / (em1 * em1);
}

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
}

bool transportMobilityDependsOnPotentials(const MobilityModelConfig& config)
{
    if (!config.jacobianFieldDerivatives)
        return false;
    return config.model == "caughey_thomas_field" ||
           config.model == "masetti_field" ||
           config.model == "caughey_thomas_field_surface" ||
           isSurfaceMobilityModel(config);
}

std::string fnv1a64(std::string_view text)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char value : text) {
        hash ^= value;
        hash *= 1099511628211ULL;
    }
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << hash;
    return output.str();
}

} // namespace

CoupledDDAssembler::CoupledDDAssembler(const DeviceMesh& mesh,
                                       const MaterialDatabase& matdb,
                                       const DopingModel& doping,
                                       double Vt,
                                       double taun,
                                       double taup)
    : CoupledDDAssembler(mesh, matdb, doping, Vt, taun, taup, {}, {})
{}

CoupledDDAssembler::CoupledDDAssembler(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    double taun,
    double taup,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges)
    : CoupledDDAssembler(mesh,
                         matdb,
                         doping,
                         Vt,
                         MobilityModelConfig{},
                         recombinationModelConfig({"srh"}, taun, taup),
                         BandgapNarrowingConfig{},
                         ImpactIonizationModelConfig{},
                         std::move(fixedCharges),
                         std::move(sheetCharges))
{}

CoupledDDAssembler::CoupledDDAssembler(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    double Vt,
    const MobilityModelConfig& mobilityConfig,
    const RecombinationModelConfig& recombinationConfig,
    const BandgapNarrowingConfig& bandgapNarrowingConfig,
    const ImpactIonizationModelConfig& impactIonizationConfig,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    DDScalingSpec scaling,
    CarrierDiagonalFloorRegularizationConfig carrierDiagonalFloor,
    CarrierStatisticsConfig carrierStatistics)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , Vt_(Vt)
    , mobilityConfig_(mobilityConfig)
    , mobility_(makeMobilityModel(mobilityConfig))
    , recombinationConfig_(recombinationConfig)
    , recombination_(recombinationConfig)
    , bandgapNarrowingConfig_(bandgapNarrowingConfig)
    , impactIonizationConfig_(impactIonizationConfig)
    , impactIonization_(makeImpactIonizationModel(impactIonizationConfig))
    , impactIonizationEnabled_(impactIonizationConfig.model != "none")
    , impactIonizationCoupled_(
          impactIonizationConfig.model != "none" &&
          impactIonizationConfig.couplingMode == "self_consistent")
    , bgnEnabled_(bandgapNarrowingConfig.model != "none")
    , carrierStatistics_(std::move(carrierStatistics))
    , ni_(detail::buildValidatedEffectiveNodeNi(
          "CoupledDDAssembler",
          mesh,
          matdb,
          doping,
          bandgapNarrowingConfig,
          Vt))
    , Nc_(detail::buildNodeDensityOfStates(
          mesh, matdb, Vt * constants::q / constants::kb, true))
    , Nv_(detail::buildNodeDensityOfStates(
          mesh, matdb, Vt * constants::q / constants::kb, false))
    , cellMaterials_(detail::buildCellMaterials(
          mesh,
          matdb,
          Vt * constants::q / constants::kb))
    , edgeCells_(detail::buildEdgeCellMap(mesh))
    , nodeCells_(detail::buildNodeCellMap(mesh))
    , vol_(detail::computeNodeVolumes(mesh))
    , couple_(detail::computeEdgeCouplings(mesh))
    , fixedInterfaceChargeRhs_(detail::computeFixedAndInterfaceChargeRhs(
          mesh,
          edgeCells_,
          fixedCharges,
          sheetCharges,
          "CoupledDDAssembler",
          scaling.enabled ? scaling.chargeAreaFactor : 1.0,
          scaling.enabled ? scaling.chargeLineFactor : 1.0))
    , scaling_(scaling)
    , carrierDiagonalFloor_(carrierDiagonalFloor)
{
    (void)usesFermiDirac(carrierStatistics_);
    if (usesFermiDirac(carrierStatistics_)) {
        for (Index node = 0; node < mesh_.numNodes(); ++node) {
            if (ni_[node] > 0.0 && (!(Nc_[node] > 0.0) || !(Nv_[node] > 0.0))) {
                throw std::invalid_argument(
                    "CoupledDDAssembler: Fermi-Dirac statistics require Nc_m3 and Nv_m3 "
                    "for every transport node.");
            }
        }
    }
    if (scaling_.enabled) {
        const auto isPositiveFinite = [](Real value) {
            return value > 0.0 && std::isfinite(value);
        };
        if (!isPositiveFinite(scaling_.V0) ||
            !isPositiveFinite(scaling_.C0) ||
            !isPositiveFinite(scaling_.mu0) ||
            !isPositiveFinite(scaling_.D0) ||
            !isPositiveFinite(scaling_.L0) ||
            !isPositiveFinite(scaling_.permittivityReference_F_per_m)) {
            throw std::invalid_argument(
                "CoupledDDAssembler: scaling references must be positive and finite when scaling is enabled.");
        }
    }
    if (carrierDiagonalFloor_.enabled) {
        if (carrierDiagonalFloor_.scale < 0.0 ||
            !std::isfinite(carrierDiagonalFloor_.scale)) {
            throw std::invalid_argument(
                "CoupledDDAssembler: carrier diagonal floor scale must be non-negative and finite.");
        }
        if (carrierDiagonalFloor_.minorityDensityRatio < 0.0 ||
            !std::isfinite(carrierDiagonalFloor_.minorityDensityRatio)) {
            throw std::invalid_argument(
                "CoupledDDAssembler: carrier diagonal floor minority density ratio must be non-negative and finite.");
        }
    }
}

VectorXd CoupledDDAssembler::pack(const CoupledDDState& state) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (state.psi.size() != N || state.phin.size() != N || state.phip.size() != N)
        throw std::invalid_argument("CoupledDDAssembler::pack: state vector size mismatch.");

    VectorXd x(3 * N);
    x.segment(psiOffset(), N) = state.psi;
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    x.segment(phinOffset(), N) =
        state.phin.array() - electronQfReference_V_ / potentialScale;
    x.segment(phipOffset(), N) =
        state.phip.array() - holeQfReference_V_ / potentialScale;
    return x;
}

CoupledDDState CoupledDDAssembler::unpack(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::unpack: vector size mismatch.");

    CoupledDDState state;
    state.psi = x.segment(psiOffset(), N);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    state.phin = x.segment(phinOffset(), N).array()
        + electronQfReference_V_ / potentialScale;
    state.phip = x.segment(phipOffset(), N).array()
        + holeQfReference_V_ / potentialScale;
    return state;
}

void CoupledDDAssembler::setQuasiFermiReferences(
    Real electronReference_V,
    Real holeReference_V)
{
    if (!std::isfinite(electronReference_V) ||
        !std::isfinite(holeReference_V)) {
        throw std::invalid_argument(
            "CoupledDDAssembler::setQuasiFermiReferences: references must be finite.");
    }
    electronQfReference_V_ = electronReference_V;
    holeQfReference_V_ = holeReference_V;
}

VectorXd CoupledDDAssembler::electronDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd n(N);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    for (int i = 0; i < N; ++i) {
        const Real psiRelative =
            x(psiOffset() + i) * potentialScale - electronQfReference_V_;
        const Real phinIncrement = x(phinOffset() + i) * potentialScale;
        const Index node = static_cast<Index>(i);
        n(i) = vela::electronDensity(
            ni_[node], Nc_[node], psiRelative, phinIncrement, Vt_,
            carrierStatistics_);
    }
    return n;
}

VectorXd CoupledDDAssembler::holeDensity(const VectorXd& x) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    VectorXd p(N);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    for (int i = 0; i < N; ++i) {
        const Real psiRelative =
            x(psiOffset() + i) * potentialScale - holeQfReference_V_;
        const Real phipIncrement = x(phipOffset() + i) * potentialScale;
        const Index node = static_cast<Index>(i);
        p(i) = vela::holeDensity(
            ni_[node], Nv_[node], psiRelative, phipIncrement, Vt_,
            carrierStatistics_);
    }
    return p;
}

bool CoupledDDAssembler::hasPositiveFiniteCarriers(const VectorXd& x) const
{
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    for (int i = 0; i < n.size(); ++i) {
        if (!std::isfinite(n(i)) || !std::isfinite(p(i)))
            return false;
        // Pure insulator nodes have ni = 0 and therefore exactly zero carrier
        // density.  Their quasi-Fermi rows are pinned by residualImpl(), so
        // zero is a valid state there.  Semiconductor and shared interface
        // nodes retain ni > 0 and must keep strictly positive carriers.
        const bool transportCarrierNode = ni_[static_cast<Index>(i)] > 0.0;
        if (transportCarrierNode && (n(i) <= 0.0 || p(i) <= 0.0))
            return false;
    }
    return true;
}

VectorXd CoupledDDAssembler::residual(const VectorXd& x,
                                      const CoupledDDBoundaryConditions& bcs) const
{
    return residualImpl(x, bcs, nullptr);
}

VectorXd CoupledDDAssembler::feedbackSubstitutionResidual(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    const CoupledDDFeedbackStateSubstitution& substitution) const
{
    return residualImpl(x, bcs, &substitution);
}

VectorXd CoupledDDAssembler::residualImpl(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    const CoupledDDFeedbackStateSubstitution* substitution) const
{
    ScopedPerformanceTimer timer("dd.residual");
    incrementPerformanceCounter("dd.residual_calls");
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::residualImpl: vector size mismatch.");

    // n and p are needed for Poisson source and configured recombination.
    VectorXd n = electronDensity(x);
    VectorXd p = holeDensity(x);
    const auto validateSubstitutionVector = [N](
        const VectorXd& values,
        const char* name,
        bool requirePositive) {
        if (values.size() != N) {
            throw std::invalid_argument(
                std::string("CoupledDDAssembler feedback substitution ") + name +
                " size mismatch.");
        }
        for (int i = 0; i < N; ++i) {
            if (!std::isfinite(values(i)) ||
                (requirePositive && values(i) <= 0.0)) {
                throw std::invalid_argument(
                    std::string("CoupledDDAssembler feedback substitution ") + name +
                    " must be finite" + (requirePositive ? " and positive." : "."));
            }
        }
    };
    if (substitution != nullptr && substitution->replaceElectronDensity) {
        validateSubstitutionVector(
            substitution->electronDensity, "electron density", true);
        n = substitution->electronDensity;
    }
    if (substitution != nullptr && substitution->replaceHoleDensity) {
        validateSubstitutionVector(
            substitution->holeDensity, "hole density", true);
        p = substitution->holeDensity;
    }
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;

    const std::vector<Real> nodeElectricFields =
        impactIonizationCoupled_
        ? detail::computeNodeElectricFields(psi, mesh_, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> bandToBandSourceIntegrals =
        recombination_.bandToBandEnabled()
        ? detail::bandToBandGenerationNodeSourceIntegrals(
              recombination_.bandToBand(), scaling_.unitSystem, mesh_,
              cellMaterials_, psi, fieldFactor)
        : std::vector<Real>{};
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig_);
    VectorXd phinIncrement = x.segment(phinOffset(), N) * potentialScale;
    VectorXd phipIncrement = x.segment(phipOffset(), N) * potentialScale;
    VectorXd phinPhysical =
        phinIncrement.array() + electronQfReference_V_;
    VectorXd phipPhysical =
        phipIncrement.array() + holeQfReference_V_;
    if (substitution != nullptr && substitution->replaceElectronQuasiFermi) {
        validateSubstitutionVector(
            substitution->electronQuasiFermi_V,
            "electron quasi-Fermi potential",
            false);
        phinPhysical = substitution->electronQuasiFermi_V;
        phinIncrement =
            phinPhysical.array() - electronQfReference_V_;
    }
    if (substitution != nullptr && substitution->replaceHoleQuasiFermi) {
        validateSubstitutionVector(
            substitution->holeQuasiFermi_V,
            "hole quasi-Fermi potential",
            false);
        phipPhysical = substitution->holeQuasiFermi_V;
        phipIncrement =
            phipPhysical.array() - holeQfReference_V_;
    }
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phinPhysical, n, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phipPhysical, p, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phinPhysical, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phipPhysical, fieldFactor)
        : std::vector<Real>{};
    const Real chargeAreaFactor = scaling_.enabled ? scaling_.chargeAreaFactor : 1.0;
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor()
        : 1.0;
    const bool sgCurrentAvalanche = impactIonizationCoupled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const std::vector<Real> sgAvalancheSourceIntegrals = sgCurrentAvalanche
        ? detail::currentDensityAvalancheSourceIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials_,
            psi,
            phinPhysical,
            phipPhysical,
            n,
            p,
            ni_,
            Vt_,
            fieldFactor,
            carrierStatistics_,
            Nc_,
            Nv_)
        : std::vector<Real>{};

    VectorXd r = VectorXd::Zero(3 * N);
    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    const auto& edgeCells = edgeCells_;
    const auto& vol = vol_;
    const auto& couple = couple_;

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = phinIncrement(i);
        const Real phin_j = phinIncrement(j);
        const Real phip_i = phipIncrement(i);
        const Real phip_j = phipIncrement(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real electronMobilityField =
            vectorQfMobility ? electronVectorMobilityFields[e] :
            (qfMobility ? std::abs((phin_j - phin_i) / h) * fieldFactor : electricField);
        const Real holeMobilityField =
            vectorQfMobility ? holeVectorMobilityFields[e] :
            (qfMobility ? std::abs((phip_j - phip_i) / h) * fieldFactor : electricField);

        const Real eps = detail::edgeEpsilon(edgeCells, mesh_, matdb_, e);
        const Real G = eps * couple[e] / h;
        const Real psiFlux = G * (psi_i - psi_j);
        r(psiOffset() + i) += psiFlux;
        r(psiOffset() + j) -= psiFlux;

        const Real mun = detail::edgeMobility(
            edgeCells, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mun * Vt_ * fieldFactor * couple[e] / h;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            // With bandgap narrowing the per-node intrinsic density varies, so
            // the ni-gradient drift term must be retained (VariableNi form).
            // Without BGN, reproduce the baseline discretization: the balanced
            // quasi-Fermi form on equal-ni edges (overflow-safe, it never
            // materializes the clamped carrier densities) and the density-based
            // SG flux on material-interface edges where ni differs.
            Real nFlux;
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (psi_i - electronQfReference_V_ - phin_i) / Vt_
                    + std::log(ni_[idxI] / Nc_[idxI]);
                const Real etaJ = (psi_j - electronQfReference_V_ - phin_j) / Vt_
                    + std::log(ni_[idxJ] / Nc_[idxJ]);
                const Real driftPotential = dpsi + Vt_ * std::log(
                    (ni_[idxJ] / Nc_[idxJ]) / (ni_[idxI] / Nc_[idxI]));
                nFlux = sgElectronFermiDiracContinuityFlux(
                    n(i), n(j), etaI, etaJ, driftPotential,
                    phin_i, phin_j, Vt_, coef);
            } else if (bgnEnabled_) {
                nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - electronQfReference_V_,
                    psi_j - electronQfReference_V_,
                    phin_i, phin_j, Vt_, coef);
            } else if (ni_[idxI] == ni_[idxJ]) {
                nFlux = sgElectronContinuityFluxFromQuasiFermiStable(
                    ni_[idxI],
                    psi_i - electronQfReference_V_,
                    psi_j - electronQfReference_V_,
                    phin_i, phin_j, Vt_, coef);
            } else {
                nFlux = sgElectronContinuityFlux(n(i), n(j), dpsi, Vt_, coef);
            }
            r(phinOffset() + i) += nFlux;
            r(phinOffset() + j) -= nFlux;
        }

        const Real mup = detail::edgeMobility(
            edgeCells, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            const Real coef = mup * Vt_ * fieldFactor * couple[e] / h;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            // See the electron flux above for the BGN gating rationale.
            Real pFlux;
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (phip_i - (psi_i - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxI] / Nv_[idxI]);
                const Real etaJ = (phip_j - (psi_j - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxJ] / Nv_[idxJ]);
                const Real driftPotential = dpsi + Vt_ * std::log(
                    (ni_[idxI] / Nv_[idxI]) / (ni_[idxJ] / Nv_[idxJ]));
                pFlux = sgHoleFermiDiracContinuityFlux(
                    p(i), p(j), etaI, etaJ, driftPotential,
                    phip_i, phip_j, Vt_, coef);
            } else if (bgnEnabled_) {
                pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - holeQfReference_V_,
                    psi_j - holeQfReference_V_,
                    phip_i, phip_j, Vt_, coef);
            } else if (ni_[idxI] == ni_[idxJ]) {
                pFlux = sgHoleContinuityFluxFromQuasiFermiStable(
                    ni_[idxI],
                    psi_i - holeQfReference_V_,
                    psi_j - holeQfReference_V_,
                    phip_i, phip_j, Vt_, coef);
            } else {
                pFlux = sgHoleContinuityFlux(p(i), p(j), dpsi, Vt_, coef);
            }
            r(phipOffset() + i) += pFlux;
            r(phipOffset() + j) -= pFlux;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        r(psiOffset() + ii) -= constants::q *
            (p(ii) - n(ii) + doping_.netDoping(i)) * vol[i] * chargeAreaFactor;

        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            // Compute n*p - ni^2 via the identity n*p = ni^2 * exp((phip-phin)/Vt),
            // i.e. n*p - ni^2 = ni^2 * expm1((phip-phin)/Vt). This avoids
            // catastrophic cancellation when phip ~= phin (near equilibrium), where
            // the naive form n*p - ni^2 is dominated by floating-point rounding
            // in exp(+u) * exp(-u) != 1.
            const bool substitutedQuasiFermi =
                substitution != nullptr &&
                (substitution->replaceElectronQuasiFermi ||
                 substitution->replaceHoleQuasiFermi);
            const Real dPhi = substitutedQuasiFermi
                ? phipPhysical(ii) - phinPhysical(ii)
                : (holeQfReference_V_ - electronQfReference_V_)
                    + (x(phipOffset() + ii) - x(phinOffset() + ii))
                        * potentialScale;
            Real recombinationNi = ni;
            Real excessProduct = 0.0;
            if (usesFermiDirac(carrierStatistics_)) {
                if (dPhi != 0.0) {
                    const Real equilibriumProduct = equilibriumCarrierProduct(
                        n(ii), p(ii), ni, Nc_[i], Nv_[i], Vt_, carrierStatistics_);
                    recombinationNi = std::sqrt(std::max<Real>(equilibriumProduct, 0.0));
                    excessProduct = n(ii) * p(ii) - equilibriumProduct;
                }
            } else {
                excessProduct = ni * ni * std::expm1(dPhi / Vt_);
            }
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), recombinationNi);
            if (R != 0.0) {
                r(phinOffset() + ii) += R * vol[i] * sourceIntegralFactor;
                r(phipOffset() + ii) += R * vol[i] * sourceIntegralFactor;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (recombination_.bandToBandEnabled()) {
            const Real source = bandToBandSourceIntegrals[i] * sourceIntegralFactor;
            if (source != 0.0) {
                r(phinOffset() + ii) -= source;
                r(phipOffset() + ii) -= source;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (impactIonizationCoupled_ && sgCurrentAvalanche) {
            const Real source = sgAvalancheSourceIntegrals[i] * sourceIntegralFactor;
            if (source != 0.0) {
                r(phinOffset() + ii) -= source;
                r(phipOffset() + ii) -= source;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        } else if (impactIonizationCoupled_) {
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                r(phinOffset() + ii) -= G * vol[i] * sourceIntegralFactor;
                r(phipOffset() + ii) -= G * vol[i] * sourceIntegralFactor;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    // Insulating nodes such as SiO2 (mun = mup = ni = 0) can have no
    // transport or recombination contribution in the continuity equations.
    // Pin those otherwise unconstrained quasi-Fermi unknowns to avoid zero
    // residual/Jacobian rows. Explicit boundary conditions below take
    // precedence over this internal gauge constraint.
    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            r(phinOffset() + ii) = x(phinOffset() + ii);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            r(phipOffset() + ii) = x(phipOffset() + ii);
    }

    for (int i = 0; i < N; ++i)
        r(psiOffset() + i) -= fixedInterfaceChargeRhs_(i);

    if (scaling_.enabled) {
        const Real poissonScale =
            scaling_.permittivityReference_F_per_m * scaling_.V0;
        const Real continuityScale = scaling_.C0 * scaling_.D0;
        for (int i = 0; i < N; ++i) {
            r(psiOffset() + i) /= poissonScale;
            r(phinOffset() + i) /= continuityScale;
            r(phipOffset() + i) /= continuityScale;
        }
    }

    if (scaling_.enabled) {
        for (Index i = 0; i < Nidx; ++i) {
            const int ii = static_cast<int>(i);
            if (!hasElectronContribution[static_cast<std::size_t>(ii)])
                r(phinOffset() + ii) = x(phinOffset() + ii);
            if (!hasHoleContribution[static_cast<std::size_t>(ii)])
                r(phipOffset() + ii) = x(phipOffset() + ii);
        }
    }

    // Boundary-condition maps are independent so multi-terminal MOS callers can
    // pin electrostatic potential, electron quasi-Fermi potential, and hole
    // quasi-Fermi potential on the source/drain/body nodes selected by contact
    // name without relying on contact order.
    for (const auto& [node, value] : bcs.psi)
        r(psiOffset() + static_cast<int>(node)) = x(psiOffset() + static_cast<int>(node)) - value;
    for (const auto& [node, value] : bcs.phin)
        r(phinOffset() + static_cast<int>(node)) =
            x(phinOffset() + static_cast<int>(node))
            - (value - electronQfReference_V_ / potentialScale);
    for (const auto& [node, value] : bcs.phip)
        r(phipOffset() + static_cast<int>(node)) =
            x(phipOffset() + static_cast<int>(node))
            - (value - holeQfReference_V_ / potentialScale);

    return r;
}

std::vector<CoupledDDCarrierTermDiagnostic>
CoupledDDAssembler::carrierContinuityTermDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    return carrierContinuityTermDiagnosticsImpl(
        x, bcs, nullptr, impactIonizationEnabled_);
}

std::vector<CoupledDDCarrierTermDiagnostic>
CoupledDDAssembler::carrierContinuityEquationTermDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    return carrierContinuityTermDiagnosticsImpl(
        x, bcs, nullptr, impactIonizationCoupled_);
}

std::vector<CoupledDDCarrierTermDiagnostic>
CoupledDDAssembler::feedbackSubstitutionCarrierContinuityTermDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    const CoupledDDFeedbackStateSubstitution& substitution) const
{
    return carrierContinuityTermDiagnosticsImpl(
        x, bcs, &substitution, impactIonizationEnabled_);
}

std::vector<CoupledDDCarrierTermDiagnostic>
CoupledDDAssembler::carrierContinuityTermDiagnosticsImpl(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    const CoupledDDFeedbackStateSubstitution* substitution,
    bool includeImpactIonization) const
{
    ScopedPerformanceTimer timer("dd.continuity_diagnostics");
    incrementPerformanceCounter("dd.continuity_diagnostics_calls");
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument(
            "CoupledDDAssembler::carrierContinuityTermDiagnostics: vector size mismatch.");

    std::vector<CoupledDDCarrierTermDiagnostic> terms(Nidx);
    for (Index i = 0; i < Nidx; ++i)
        terms[i].nodeId = i;

    VectorXd n = electronDensity(x);
    VectorXd p = holeDensity(x);
    const auto validateSubstitutionVector = [N](
        const VectorXd& values,
        const char* name,
        bool requirePositive) {
        if (values.size() != N) {
            throw std::invalid_argument(
                std::string("CoupledDDAssembler feedback substitution ") + name +
                " size mismatch.");
        }
        for (int i = 0; i < N; ++i) {
            if (!std::isfinite(values(i)) ||
                (requirePositive && values(i) <= 0.0)) {
                throw std::invalid_argument(
                    std::string("CoupledDDAssembler feedback substitution ") + name +
                    " must be finite" + (requirePositive ? " and positive." : "."));
            }
        }
    };
    if (substitution != nullptr && substitution->replaceElectronDensity) {
        validateSubstitutionVector(
            substitution->electronDensity, "electron density", true);
        n = substitution->electronDensity;
    }
    if (substitution != nullptr && substitution->replaceHoleDensity) {
        validateSubstitutionVector(
            substitution->holeDensity, "hole density", true);
        p = substitution->holeDensity;
    }
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;

    const std::vector<Real> nodeElectricFields =
        includeImpactIonization
        ? detail::computeNodeElectricFields(psi, mesh_, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> bandToBandSourceIntegrals =
        recombination_.bandToBandEnabled()
        ? detail::bandToBandGenerationNodeSourceIntegrals(
              recombination_.bandToBand(), scaling_.unitSystem, mesh_,
              cellMaterials_, psi, fieldFactor)
        : std::vector<Real>{};
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig_);
    VectorXd phinIncrement = x.segment(phinOffset(), N) * potentialScale;
    VectorXd phipIncrement = x.segment(phipOffset(), N) * potentialScale;
    VectorXd phinPhysical =
        phinIncrement.array() + electronQfReference_V_;
    VectorXd phipPhysical =
        phipIncrement.array() + holeQfReference_V_;
    if (substitution != nullptr && substitution->replaceElectronQuasiFermi) {
        validateSubstitutionVector(
            substitution->electronQuasiFermi_V,
            "electron quasi-Fermi potential",
            false);
        phinPhysical = substitution->electronQuasiFermi_V;
        phinIncrement =
            phinPhysical.array() - electronQfReference_V_;
    }
    if (substitution != nullptr && substitution->replaceHoleQuasiFermi) {
        validateSubstitutionVector(
            substitution->holeQuasiFermi_V,
            "hole quasi-Fermi potential",
            false);
        phipPhysical = substitution->holeQuasiFermi_V;
        phipIncrement =
            phipPhysical.array() - holeQfReference_V_;
    }
    const std::vector<Real> nodeElectronDrivingFields = (includeImpactIonization && qfImpact)
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phinPhysical, n, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (includeImpactIonization && qfImpact)
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phipPhysical, p, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phinPhysical, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phipPhysical, fieldFactor)
        : std::vector<Real>{};
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor()
        : 1.0;
    const bool sgCurrentAvalanche = includeImpactIonization &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const detail::SgAvalancheSourceComponentIntegrals sgAvalancheSourceComponents =
        sgCurrentAvalanche
        ? detail::currentDensityAvalancheSourceComponentIntegrals(
            impactIonizationConfig_,
            *impactIonization_,
            mobilityConfig_,
            *mobility_,
            edgeCells_,
            mesh_,
            doping_,
            cellMaterials_,
            psi,
            phinPhysical,
            phipPhysical,
            n,
            p,
            ni_,
            Vt_,
            fieldFactor,
            carrierStatistics_,
            Nc_,
            Nv_)
        : detail::SgAvalancheSourceComponentIntegrals{};

    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = phinIncrement(i);
        const Real phin_j = phinIncrement(j);
        const Real phip_i = phipIncrement(i);
        const Real phip_j = phipIncrement(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real electronMobilityField =
            vectorQfMobility ? electronVectorMobilityFields[e] :
            (qfMobility ? std::abs((phin_j - phin_i) / h) * fieldFactor : electricField);
        const Real holeMobilityField =
            vectorQfMobility ? holeVectorMobilityFields[e] :
            (qfMobility ? std::abs((phip_j - phip_i) / h) * fieldFactor : electricField);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;
            const Real coef = mun * Vt_ * fieldFactor * couple_[e] / h;
            Real nFlux = 0.0;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (psi_i - electronQfReference_V_ - phin_i) / Vt_
                    + std::log(ni_[idxI] / Nc_[idxI]);
                const Real etaJ = (psi_j - electronQfReference_V_ - phin_j) / Vt_
                    + std::log(ni_[idxJ] / Nc_[idxJ]);
                const Real drift = dpsi + Vt_ * std::log(
                    (ni_[idxJ] / Nc_[idxJ]) / (ni_[idxI] / Nc_[idxI]));
                nFlux = sgElectronFermiDiracContinuityFlux(
                    n(i), n(j), etaI, etaJ, drift, phin_i, phin_j, Vt_, coef);
            } else {
                nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - electronQfReference_V_,
                    psi_j - electronQfReference_V_,
                    phin_i, phin_j, Vt_, coef, bgnEnabled_);
            }
            terms[static_cast<Index>(i)].electronFlux += nFlux;
            terms[static_cast<Index>(j)].electronFlux -= nFlux;
            terms[static_cast<Index>(i)].electronFluxAbsSum += std::abs(nFlux);
            terms[static_cast<Index>(j)].electronFluxAbsSum += std::abs(nFlux);
        }

        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;
            const Real coef = mup * Vt_ * fieldFactor * couple_[e] / h;
            Real pFlux = 0.0;
            const Index idxI = static_cast<Index>(i);
            const Index idxJ = static_cast<Index>(j);
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (phip_i - (psi_i - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxI] / Nv_[idxI]);
                const Real etaJ = (phip_j - (psi_j - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxJ] / Nv_[idxJ]);
                const Real drift = dpsi + Vt_ * std::log(
                    (ni_[idxI] / Nv_[idxI]) / (ni_[idxJ] / Nv_[idxJ]));
                pFlux = sgHoleFermiDiracContinuityFlux(
                    p(i), p(j), etaI, etaJ, drift, phip_i, phip_j, Vt_, coef);
            } else {
                pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - holeQfReference_V_,
                    psi_j - holeQfReference_V_,
                    phip_i, phip_j, Vt_, coef, bgnEnabled_);
            }
            terms[static_cast<Index>(i)].holeFlux += pFlux;
            terms[static_cast<Index>(j)].holeFlux -= pFlux;
            terms[static_cast<Index>(i)].holeFluxAbsSum += std::abs(pFlux);
            terms[static_cast<Index>(j)].holeFluxAbsSum += std::abs(pFlux);
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            const bool substitutedQuasiFermi =
                substitution != nullptr &&
                (substitution->replaceElectronQuasiFermi ||
                 substitution->replaceHoleQuasiFermi);
            const Real dPhi = substitutedQuasiFermi
                ? phipPhysical(ii) - phinPhysical(ii)
                : (holeQfReference_V_ - electronQfReference_V_)
                    + (x(phipOffset() + ii) - x(phinOffset() + ii))
                        * potentialScale;
            Real recombinationNi = ni;
            Real excessProduct = 0.0;
            if (usesFermiDirac(carrierStatistics_)) {
                if (dPhi != 0.0) {
                    const Real equilibriumProduct = equilibriumCarrierProduct(
                        n(ii), p(ii), ni, Nc_[i], Nv_[i], Vt_, carrierStatistics_);
                    recombinationNi = std::sqrt(std::max<Real>(equilibriumProduct, 0.0));
                    excessProduct = n(ii) * p(ii) - equilibriumProduct;
                }
            } else {
                excessProduct = ni * ni * std::expm1(dPhi / Vt_);
            }
            const Real R = recombination_.totalRateFromExcessProduct(
                excessProduct, n(ii), p(ii), recombinationNi);
            if (R != 0.0) {
                const Real contribution = R * vol_[i] * sourceIntegralFactor;
                terms[i].electronRecombination += contribution;
                terms[i].holeRecombination += contribution;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (recombination_.bandToBandEnabled()) {
            const Real contribution =
                -bandToBandSourceIntegrals[i] * sourceIntegralFactor;
            if (contribution != 0.0) {
                terms[i].electronRecombination += contribution;
                terms[i].holeRecombination += contribution;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }

        if (includeImpactIonization && sgCurrentAvalanche) {
            const Real source =
                sgAvalancheSourceComponents.combined[i] * sourceIntegralFactor;
            const Real contribution = -source;
            if (contribution != 0.0) {
                terms[i].electronImpact += contribution;
                terms[i].holeImpact += contribution;
                terms[i].impactElectronSource +=
                    sgAvalancheSourceComponents.electron[i] * sourceIntegralFactor;
                terms[i].impactHoleSource +=
                    sgAvalancheSourceComponents.hole[i] * sourceIntegralFactor;
                terms[i].impactCombinedSource += source;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        } else if (includeImpactIonization) {
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                const Real contribution = -G * vol_[i] * sourceIntegralFactor;
                terms[i].electronImpact += contribution;
                terms[i].holeImpact += contribution;
                terms[i].impactCombinedSource += G * vol_[i] * sourceIntegralFactor;
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    if (scaling_.enabled) {
        const Real continuityScale = scaling_.C0 * scaling_.D0;
        for (auto& term : terms) {
            term.electronFlux /= continuityScale;
            term.holeFlux /= continuityScale;
            term.electronFluxAbsSum /= continuityScale;
            term.holeFluxAbsSum /= continuityScale;
            term.electronRecombination /= continuityScale;
            term.holeRecombination /= continuityScale;
            term.electronImpact /= continuityScale;
            term.holeImpact /= continuityScale;
            term.impactElectronSource /= continuityScale;
            term.impactHoleSource /= continuityScale;
            term.impactCombinedSource /= continuityScale;
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            terms[i].electronGauge = x(phinOffset() + ii);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            terms[i].holeGauge = x(phipOffset() + ii);
    }

    for (auto& term : terms) {
        term.electronResidual = term.electronFlux
            + term.electronRecombination
            + term.electronImpact
            + term.electronGauge;
        term.holeResidual = term.holeFlux
            + term.holeRecombination
            + term.holeImpact
            + term.holeGauge;
    }

    auto applyElectronBoundary = [&](Index node, Real value) {
        CoupledDDCarrierTermDiagnostic& term = terms[node];
        term.electronFlux = 0.0;
        term.electronFluxAbsSum = 0.0;
        term.electronRecombination = 0.0;
        term.electronImpact = 0.0;
        term.impactElectronSource = 0.0;
        term.impactHoleSource = 0.0;
        term.impactCombinedSource = 0.0;
        term.electronGauge = 0.0;
        term.electronBoundary =
            x(phinOffset() + static_cast<int>(node))
            - (value - electronQfReference_V_ / potentialScale);
        term.electronResidual = term.electronBoundary;
    };
    auto applyHoleBoundary = [&](Index node, Real value) {
        CoupledDDCarrierTermDiagnostic& term = terms[node];
        term.holeFlux = 0.0;
        term.holeFluxAbsSum = 0.0;
        term.holeRecombination = 0.0;
        term.holeImpact = 0.0;
        term.impactElectronSource = 0.0;
        term.impactHoleSource = 0.0;
        term.impactCombinedSource = 0.0;
        term.holeGauge = 0.0;
        term.holeBoundary =
            x(phipOffset() + static_cast<int>(node))
            - (value - holeQfReference_V_ / potentialScale);
        term.holeResidual = term.holeBoundary;
    };
    for (const auto& [node, value] : bcs.phin)
        applyElectronBoundary(node, value);
    for (const auto& [node, value] : bcs.phip)
        applyHoleBoundary(node, value);

    return terms;
}

std::vector<CoupledDDEdgeFluxDiagnostic>
CoupledDDAssembler::sgEdgeFluxDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    (void)bcs;
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument(
            "CoupledDDAssembler::sgEdgeFluxDiagnostics: vector size mismatch.");

    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);

    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const VectorXd phinField = x.segment(phinOffset(), N) * potentialScale;
    const VectorXd phipField = x.segment(phipOffset(), N) * potentialScale;
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phinField, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phipField, fieldFactor)
        : std::vector<Real>{};
    const Real continuityScale =
        scaling_.enabled ? (scaling_.C0 * scaling_.D0) : 1.0;

    std::vector<CoupledDDEdgeFluxDiagnostic> edges;
    edges.reserve(static_cast<std::size_t>(mesh_.numEdges()));

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real electronMobilityField =
            vectorQfMobility ? electronVectorMobilityFields[e] :
            (qfMobility ? std::abs((phin_j - phin_i) / h) * fieldFactor : electricField);
        const Real holeMobilityField =
            vectorQfMobility ? holeVectorMobilityFields[e] :
            (qfMobility ? std::abs((phip_j - phip_i) / h) * fieldFactor : electricField);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, &psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, &psi);

        Real nFlux = 0.0;
        if (mun > 0.0) {
            const Real coef = mun * Vt_ * fieldFactor * couple_[e] / h;
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (psi_i - electronQfReference_V_ - phin_i) / Vt_
                    + std::log(ni_[idxI] / Nc_[idxI]);
                const Real etaJ = (psi_j - electronQfReference_V_ - phin_j) / Vt_
                    + std::log(ni_[idxJ] / Nc_[idxJ]);
                const Real drift = dpsi + Vt_ * std::log(
                    (ni_[idxJ] / Nc_[idxJ]) / (ni_[idxI] / Nc_[idxI]));
                nFlux = sgElectronFermiDiracContinuityFlux(
                    n(i), n(j), etaI, etaJ, drift,
                    phin_i, phin_j, Vt_, coef);
            } else {
                nFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - electronQfReference_V_,
                    psi_j - electronQfReference_V_,
                    phin_i, phin_j, Vt_, coef, bgnEnabled_);
            }
        }
        Real pFlux = 0.0;
        if (mup > 0.0) {
            const Real coef = mup * Vt_ * fieldFactor * couple_[e] / h;
            if (usesFermiDirac(carrierStatistics_)) {
                const Real etaI = (phip_i - (psi_i - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxI] / Nv_[idxI]);
                const Real etaJ = (phip_j - (psi_j - holeQfReference_V_)) / Vt_
                    + std::log(ni_[idxJ] / Nv_[idxJ]);
                const Real drift = dpsi + Vt_ * std::log(
                    (ni_[idxI] / Nv_[idxI]) / (ni_[idxJ] / Nv_[idxJ]));
                pFlux = sgHoleFermiDiracContinuityFlux(
                    p(i), p(j), etaI, etaJ, drift,
                    phip_i, phip_j, Vt_, coef);
            } else {
                pFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[idxI], ni_[idxJ],
                    psi_i - holeQfReference_V_,
                    psi_j - holeQfReference_V_,
                    phip_i, phip_j, Vt_, coef, bgnEnabled_);
            }
        }

        const Node& node0 = mesh_.getNode(edge.n0);
        const Node& node1 = mesh_.getNode(edge.n1);
        CoupledDDEdgeFluxDiagnostic record;
        record.edgeId = e;
        record.node0 = idxI;
        record.node1 = idxJ;
        record.x0 = node0.x;
        record.y0 = node0.y;
        record.x1 = node1.x;
        record.y1 = node1.y;
        record.length_m = h;
        record.couple_m = couple_[e];
        record.netDopingAvg_m3 =
            0.5 * (doping_.netDoping(idxI) + doping_.netDoping(idxJ));
        record.ni0_m3 = ni_[idxI];
        record.ni1_m3 = ni_[idxJ];
        record.psi0_V = psi_i;
        record.psi1_V = psi_j;
        record.phin0_V = phin_i + electronQfReference_V_;
        record.phin1_V = phin_j + electronQfReference_V_;
        record.phip0_V = phip_i + holeQfReference_V_;
        record.phip1_V = phip_j + holeQfReference_V_;
        record.electricField_V_m = electricField;
        record.electronMobility_m2_V_s = mun;
        record.holeMobility_m2_V_s = mup;
        record.electronFlux = nFlux / continuityScale;
        record.holeFlux = pFlux / continuityScale;
        record.electronParticleLineFlux_per_m_s = nFlux;
        record.holeParticleLineFlux_per_m_s = pFlux;
        edges.push_back(record);
    }

    return edges;
}

std::vector<CoupledDDTransportEdgeJacobianDiagnostic>
CoupledDDAssembler::transportEdgeJacobianDiagnostics(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real physicalFiniteDifferenceStep_V) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N)
        throw std::invalid_argument(
            "CoupledDDAssembler::transportEdgeJacobianDiagnostics: vector size mismatch.");
    if (!(physicalFiniteDifferenceStep_V > 0.0))
        throw std::invalid_argument(
            "CoupledDDAssembler::transportEdgeJacobianDiagnostics: finite-difference step must be positive.");

    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled
        ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const Real continuityScale = scaling_.enabled
        ? scaling_.C0 * scaling_.D0 : 1.0;
    const Real derivativeScale = scaling_.enabled
        ? scaling_.V0 / continuityScale : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const bool qfMobility =
        mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";

    auto isConstrained = [&](CarrierType carrier, Index node) {
        return carrier == CarrierType::Electron
            ? bcs.phin.find(node) != bcs.phin.end()
            : bcs.phip.find(node) != bcs.phip.end();
    };

    std::vector<CoupledDDTransportEdgeJacobianDiagnostic> records;
    records.reserve(static_cast<std::size_t>(mesh_.numEdges()) * 8);
    for (Index edgeId = 0; edgeId < mesh_.numEdges(); ++edgeId) {
        const Edge& edge = mesh_.getEdge(edgeId);
        const Real h = edge.length;
        if (h <= 1.0e-30 || couple_[edgeId] <= 0.0)
            continue;
        const Index nodes[2] = {edge.n0, edge.n1};
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi0 = psi(i);
        const Real psi1 = psi(j);
        const Real electricField = std::abs((psi1 - psi0) / h) * fieldFactor;

        for (const CarrierType carrier : {CarrierType::Electron, CarrierType::Hole}) {
            const bool electron = carrier == CarrierType::Electron;
            const int qfpOffset = electron ? phinOffset() : phipOffset();
            const Real qfpReference = electron
                ? electronQfReference_V_ : holeQfReference_V_;
            const Real qfp0 = x(qfpOffset + i) * potentialScale;
            const Real qfp1 = x(qfpOffset + j) * potentialScale;
            const Real qfpField = std::abs((qfp1 - qfp0) / h) * fieldFactor;
            const Real drive = qfMobility ? qfpField : electricField;
            const Real mobility = detail::edgeMobility(
                edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, edgeId,
                carrier, drive, &mobilityConfig_, &psi);
            if (mobility <= 0.0)
                continue;

            auto mobilityAt = [&](Real endpoint0, Real endpoint1) {
                const Real localQfpField =
                    std::abs((endpoint1 - endpoint0) / h) * fieldFactor;
                const Real localDrive = qfMobility ? localQfpField : electricField;
                return detail::edgeMobility(
                    edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, edgeId,
                    carrier, localDrive, &mobilityConfig_, &psi);
            };
            auto fluxAt = [&](Real endpoint0, Real endpoint1, Real localMobility) {
                const Real coef =
                    localMobility * Vt_ * fieldFactor * couple_[edgeId] / h;
                if (electron) {
                    return sgElectronContinuityFluxFromQuasiFermiVariableNi(
                        ni_[edge.n0], ni_[edge.n1],
                        psi0 - qfpReference, psi1 - qfpReference,
                        endpoint0, endpoint1, Vt_, coef, bgnEnabled_);
                }
                return sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni_[edge.n0], ni_[edge.n1],
                    psi0 - qfpReference, psi1 - qfpReference,
                    endpoint0, endpoint1, Vt_, coef, bgnEnabled_);
            };

            const Real coefficient =
                mobility * Vt_ * fieldFactor * couple_[edgeId] / h;
            const Real eta = electron
                ? (psi1 - psi0) / Vt_ + std::log(ni_[edge.n1] / ni_[edge.n0])
                : (psi1 - psi0) / Vt_ + std::log(ni_[edge.n0] / ni_[edge.n1]);
            const Real bPlus = bernoulli(eta);
            const Real bMinus = bernoulli(-eta);
            const Real density0 = electron
                ? ni_[edge.n0] * limitedExp((psi0 - qfpReference - qfp0) / Vt_)
                : ni_[edge.n0] * limitedExp((qfp0 - (psi0 - qfpReference)) / Vt_);
            const Real density1 = electron
                ? ni_[edge.n1] * limitedExp((psi1 - qfpReference - qfp1) / Vt_)
                : ni_[edge.n1] * limitedExp((qfp1 - (psi1 - qfpReference)) / Vt_);
            const Real flux = fluxAt(qfp0, qfp1, mobility);

            for (int columnEndpoint = 0; columnEndpoint < 2; ++columnEndpoint) {
                Real plus0 = qfp0;
                Real plus1 = qfp1;
                Real minus0 = qfp0;
                Real minus1 = qfp1;
                if (columnEndpoint == 0) {
                    plus0 += physicalFiniteDifferenceStep_V;
                    minus0 -= physicalFiniteDifferenceStep_V;
                } else {
                    plus1 += physicalFiniteDifferenceStep_V;
                    minus1 -= physicalFiniteDifferenceStep_V;
                }
                const Real mobilityPlus = mobilityAt(plus0, plus1);
                const Real mobilityMinus = mobilityAt(minus0, minus1);
                const Real liveDerivative =
                    (fluxAt(plus0, plus1, mobilityPlus)
                     - fluxAt(minus0, minus1, mobilityMinus))
                    / (2.0 * physicalFiniteDifferenceStep_V);
                const Real frozenDerivative =
                    (fluxAt(plus0, plus1, mobility)
                     - fluxAt(minus0, minus1, mobility))
                    / (2.0 * physicalFiniteDifferenceStep_V);
                const Real mobilityDerivative =
                    (mobilityPlus - mobilityMinus)
                    / (2.0 * physicalFiniteDifferenceStep_V);
                const Real drivePlus = qfMobility
                    ? std::abs((plus1 - plus0) / h) * fieldFactor : electricField;
                const Real driveMinus = qfMobility
                    ? std::abs((minus1 - minus0) / h) * fieldFactor : electricField;
                const Real driveDerivative =
                    (drivePlus - driveMinus)
                    / (2.0 * physicalFiniteDifferenceStep_V);
                const Real productionDerivative = electron
                    ? (columnEndpoint == 0
                        ? coefficient * bMinus * (-density0 / Vt_)
                        : coefficient * bPlus * (density1 / Vt_))
                    : (columnEndpoint == 0
                        ? coefficient * bPlus * (density0 / Vt_)
                        : coefficient * bMinus * (-density1 / Vt_));

                for (int rowEndpoint = 0; rowEndpoint < 2; ++rowEndpoint) {
                    const Real rowSign = rowEndpoint == 0 ? 1.0 : -1.0;
                    const bool rowConstrained = isConstrained(carrier, nodes[rowEndpoint]);
                    const bool columnConstrained =
                        isConstrained(carrier, nodes[columnEndpoint]);
                    CoupledDDTransportEdgeJacobianDiagnostic record;
                    record.edgeId = edgeId;
                    record.carrier = electron ? "electron" : "hole";
                    record.node0 = edge.n0;
                    record.node1 = edge.n1;
                    record.rowNode = nodes[rowEndpoint];
                    record.columnNode = nodes[columnEndpoint];
                    record.rowEndpoint = rowEndpoint;
                    record.columnEndpoint = columnEndpoint;
                    record.lengthInternal = h;
                    record.coupleInternal = couple_[edgeId];
                    record.qfpDriveInternal = drive;
                    record.dQfpDriveDColumnInternal = driveDerivative;
                    record.mobilityInternal = mobility;
                    record.dMobilityDColumnInternal = mobilityDerivative;
                    record.bernoulliNode0 = electron ? bMinus : bPlus;
                    record.bernoulliNode1 = electron ? bPlus : bMinus;
                    record.carrierDensityNode0Internal = density0;
                    record.carrierDensityNode1Internal = density1;
                    record.fluxPhysical = flux;
                    record.fluxScaled = flux / continuityScale;
                    record.rowSign = rowSign;
                    record.productionFrozenMobilityDerivativePhysical =
                        productionDerivative;
                    record.frozenMobilityFiniteDifferenceDerivativePhysical =
                        frozenDerivative;
                    record.liveMobilityFiniteDifferenceDerivativePhysical =
                        liveDerivative;
                    record.mobilityResponseDerivativePhysical =
                        flux / mobility * mobilityDerivative;
                    record.liveMinusFrozenFiniteDifferenceDerivativePhysical =
                        liveDerivative - frozenDerivative;
                    // The Bernoulli argument contains psi and effective-ni but
                    // no same-carrier QFP.  Its QFP-column derivative is exactly zero.
                    record.bernoulliQfpDerivativePhysical = 0.0;
                    record.carrierPopulationDerivativePhysical =
                        productionDerivative;
                    record.productionRowDerivativeScaled =
                        rowSign * productionDerivative * derivativeScale;
                    record.liveMobilityRowDerivativeScaled =
                        rowSign * liveDerivative * derivativeScale;
                    record.rowConstrained = rowConstrained;
                    record.columnConstrained = columnConstrained;
                    record.contactEliminatedProductionEdgeDerivative =
                        rowConstrained ? 0.0 : record.productionRowDerivativeScaled;
                    record.contactIdentityEntry =
                        rowConstrained && record.rowNode == record.columnNode ? 1.0 : 0.0;
                    records.push_back(std::move(record));
                }
            }
        }
    }
    return records;
}

VectorXd CoupledDDAssembler::impactIonizationSourceResidual(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceResidual: vector size mismatch.");
    }

    VectorXd source = VectorXd::Zero(3 * N);
    const auto terms = carrierContinuityTermDiagnostics(x, bcs);
    for (int node = 0; node < N; ++node) {
        source(phinOffset() + node) =
            terms[static_cast<std::size_t>(node)].electronImpact;
        source(phipOffset() + node) =
            terms[static_cast<std::size_t>(node)].holeImpact;
    }
    return source;
}

SparseMatrixd CoupledDDAssembler::impactIonizationSourceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceJacobian: vector size mismatch.");
    }
    if (!impactIonizationEnabled_ ||
        !detail::usesElementEdgeGssLauxAvalancheSource(impactIonizationConfig_)) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceJacobian requires the "
            "element_edge_sg_gss_laux avalanche source.");
    }
    if (isSurfaceMobilityModel(mobilityConfig_)) {
        throw std::invalid_argument(
            "element-edge source-only Jacobian does not support surface mobility");
    }

    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled
        ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor() : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const VectorXd phin =
        (x.segment(phinOffset(), N) * potentialScale).array()
        + electronQfReference_V_;
    const VectorXd phip =
        (x.segment(phipOffset(), N) * potentialScale).array()
        + holeQfReference_V_;
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells_, mesh_);

    std::vector<bool> constrainedRows(static_cast<std::size_t>(3 * N), false);
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(
            phinOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(
            phipOffset() + static_cast<int>(node))] = true;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(mesh_.numCells()) * 54);
    const Real derivativeScale = scaling_.enabled
        ? scaling_.V0 / (scaling_.C0 * scaling_.D0) : 1.0;

    for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
        const Cell& cell = mesh_.getCell(cellId);
        std::array<detail::Tri3LocalForwardDual, 3> localPsi{};
        std::array<detail::Tri3LocalForwardDual, 3> localPhin{};
        std::array<detail::Tri3LocalForwardDual, 3> localPhip{};
        std::array<detail::Tri3LocalForwardDual, 3> localElectronDensity{};
        std::array<detail::Tri3LocalForwardDual, 3> localHoleDensity{};
        std::array<Real, 3> localIntrinsicDensity{};
        for (std::size_t localNode = 0; localNode < 3; ++localNode) {
            const Index node = cell.node_ids[localNode];
            const int nodeIndex = static_cast<int>(node);
            localPsi[localNode] = detail::Tri3LocalForwardDual::variable(
                psi(nodeIndex), localNode);
            localPhin[localNode] = detail::Tri3LocalForwardDual::variable(
                phin(nodeIndex), 3 + localNode);
            localPhip[localNode] = detail::Tri3LocalForwardDual::variable(
                phip(nodeIndex), 6 + localNode);
            localIntrinsicDensity[localNode] = ni_[node];
            localElectronDensity[localNode] =
                detail::Tri3LocalForwardDual(ni_[node]) * detail::localAdLimitedExp(
                    (localPsi[localNode] - localPhin[localNode]) /
                    detail::Tri3LocalForwardDual(Vt_));
            localHoleDensity[localNode] =
                detail::Tri3LocalForwardDual(ni_[node]) * detail::localAdLimitedExp(
                    (localPhip[localNode] - localPsi[localNode]) /
                    detail::Tri3LocalForwardDual(Vt_));
        }
        const auto sources =
            detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<
                detail::Tri3LocalForwardDual>(
                impactIonizationConfig_, mobilityConfig_, *mobility_,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh_, doping_,
                cellMaterials_, cellId, localPsi, localPhin, localPhip,
                localElectronDensity, localHoleDensity, localIntrinsicDensity,
                Vt_, fieldFactor);

        for (int localRow = 0; localRow < 3; ++localRow) {
            const int rowNode = static_cast<int>(
                cell.node_ids[static_cast<std::size_t>(localRow)]);
            for (std::size_t localDof = 0;
                 localDof < detail::Tri3LocalPotentialDofCount; ++localDof) {
                const int variableBlock = static_cast<int>(localDof / 3);
                const int localColumn = static_cast<int>(localDof % 3);
                const int columnNode = static_cast<int>(
                    cell.node_ids[static_cast<std::size_t>(localColumn)]);
                const int columnOffset = variableBlock == 0
                    ? psiOffset()
                    : (variableBlock == 1 ? phinOffset() : phipOffset());
                const Real derivative =
                    -sources.combined[static_cast<std::size_t>(localRow)]
                         .derivative[localDof] * sourceIntegralFactor * derivativeScale;
                const int electronRow = phinOffset() + rowNode;
                const int holeRow = phipOffset() + rowNode;
                if (derivative != 0.0 &&
                    !constrainedRows[static_cast<std::size_t>(electronRow)]) {
                    triplets.emplace_back(
                        electronRow, columnOffset + columnNode, derivative);
                }
                if (derivative != 0.0 &&
                    !constrainedRows[static_cast<std::size_t>(holeRow)]) {
                    triplets.emplace_back(
                        holeRow, columnOffset + columnNode, derivative);
                }
            }
        }
    }

    SparseMatrixd jacobian(3 * N, 3 * N);
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    return jacobian;
}

SparseMatrixd CoupledDDAssembler::impactIonizationSourceFiniteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    return impactIonizationSourceFiniteDifferenceJacobianImpl<Real>(
        x, bcs, relativeStep);
}

SparseMatrixd
CoupledDDAssembler::impactIonizationSourceBranchResolvedFiniteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    using BranchResolvedScalar =
        boost::multiprecision::number<
            boost::multiprecision::cpp_dec_float<50>,
            boost::multiprecision::et_off>;
    return impactIonizationSourceFiniteDifferenceJacobianImpl<
        BranchResolvedScalar>(x, bcs, relativeStep);
}

template <typename Scalar>
SparseMatrixd
CoupledDDAssembler::impactIonizationSourceFiniteDifferenceJacobianImpl(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    const int N = static_cast<int>(mesh_.numNodes());
    if (x.size() != 3 * N) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceFiniteDifferenceJacobianImpl: "
            "vector size mismatch.");
    }
    if (relativeStep <= 0.0 || !std::isfinite(relativeStep)) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceFiniteDifferenceJacobianImpl: "
            "step must be positive and finite.");
    }
    if (!impactIonizationEnabled_ ||
        !detail::usesElementEdgeGssLauxAvalancheSource(impactIonizationConfig_)) {
        throw std::invalid_argument(
            "CoupledDDAssembler::impactIonizationSourceFiniteDifferenceJacobian "
            "requires the element_edge_sg_gss_laux avalanche source.");
    }
    if (isSurfaceMobilityModel(mobilityConfig_)) {
        throw std::invalid_argument(
            "element-edge source-only FD Jacobian does not support surface mobility");
    }

    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled
        ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor() : 1.0;
    const Real derivativeScale = scaling_.enabled
        ? scaling_.V0 / (scaling_.C0 * scaling_.D0) : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const VectorXd phin =
        (x.segment(phinOffset(), N) * potentialScale).array()
        + electronQfReference_V_;
    const VectorXd phip =
        (x.segment(phipOffset(), N) * potentialScale).array()
        + holeQfReference_V_;
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells_, mesh_);

    std::vector<bool> constrainedRows(static_cast<std::size_t>(3 * N), false);
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(
            phinOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(
            phipOffset() + static_cast<int>(node))] = true;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(mesh_.numCells()) * 54);
    for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
        const Cell& cell = mesh_.getCell(cellId);
        std::array<Scalar, 3> localPsi{};
        std::array<Scalar, 3> localPhin{};
        std::array<Scalar, 3> localPhip{};
        std::array<Real, 3> localIntrinsicDensity{};
        for (std::size_t localNode = 0; localNode < 3; ++localNode) {
            const Index node = cell.node_ids[localNode];
            const int nodeIndex = static_cast<int>(node);
            localPsi[localNode] = psi(nodeIndex);
            localPhin[localNode] = phin(nodeIndex);
            localPhip[localNode] = phip(nodeIndex);
            localIntrinsicDensity[localNode] = ni_[node];
        }
        auto evaluate = [&](const std::array<Scalar, 3>& psiValues,
                            const std::array<Scalar, 3>& phinValues,
                            const std::array<Scalar, 3>& phipValues) {
            std::array<Scalar, 3> electronDensity{};
            std::array<Scalar, 3> holeDensity{};
            for (std::size_t localNode = 0; localNode < 3; ++localNode) {
                electronDensity[localNode] =
                    Scalar(localIntrinsicDensity[localNode]) *
                    detail::localAdLimitedExp(
                        (psiValues[localNode] - phinValues[localNode]) /
                        Scalar(Vt_));
                holeDensity[localNode] =
                    Scalar(localIntrinsicDensity[localNode]) *
                    detail::localAdLimitedExp(
                        (phipValues[localNode] - psiValues[localNode]) /
                        Scalar(Vt_));
            }
            return detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<Scalar>(
                impactIonizationConfig_, mobilityConfig_, *mobility_,
                cellEdges.at(static_cast<std::size_t>(cellId)), mesh_, doping_,
                cellMaterials_, cellId, psiValues, phinValues, phipValues,
                electronDensity, holeDensity, localIntrinsicDensity, Vt_,
                fieldFactor).combined;
        };

        for (int variableBlock = 0; variableBlock < 3; ++variableBlock) {
            for (int localColumn = 0; localColumn < 3; ++localColumn) {
                auto psiPlus = localPsi;
                auto psiMinus = localPsi;
                auto phinPlus = localPhin;
                auto phinMinus = localPhin;
                auto phipPlus = localPhip;
                auto phipMinus = localPhip;
                const Scalar value = variableBlock == 0
                    ? localPsi[static_cast<std::size_t>(localColumn)]
                    : (variableBlock == 1
                        ? localPhin[static_cast<std::size_t>(localColumn)]
                        : localPhip[static_cast<std::size_t>(localColumn)]);
                const Real physicalValue = static_cast<Real>(value);
                const Scalar step = Scalar(
                    detail::physicalPotentialCentralDifferenceStep(
                        physicalValue, potentialScale, relativeStep));
                if (variableBlock == 0) {
                    psiPlus[static_cast<std::size_t>(localColumn)] += step;
                    psiMinus[static_cast<std::size_t>(localColumn)] -= step;
                } else if (variableBlock == 1) {
                    phinPlus[static_cast<std::size_t>(localColumn)] += step;
                    phinMinus[static_cast<std::size_t>(localColumn)] -= step;
                } else {
                    phipPlus[static_cast<std::size_t>(localColumn)] += step;
                    phipMinus[static_cast<std::size_t>(localColumn)] -= step;
                }
                const auto plus = evaluate(psiPlus, phinPlus, phipPlus);
                const auto minus = evaluate(psiMinus, phinMinus, phipMinus);
                const int columnNode = static_cast<int>(
                    cell.node_ids[static_cast<std::size_t>(localColumn)]);
                const int columnOffset = variableBlock == 0
                    ? psiOffset()
                    : (variableBlock == 1 ? phinOffset() : phipOffset());
                for (int localRow = 0; localRow < 3; ++localRow) {
                    const Real derivative =
                        static_cast<Real>(
                            -(plus[static_cast<std::size_t>(localRow)] -
                              minus[static_cast<std::size_t>(localRow)]) /
                            (Scalar(2.0) * step) *
                            Scalar(sourceIntegralFactor * derivativeScale));
                    const int rowNode = static_cast<int>(
                        cell.node_ids[static_cast<std::size_t>(localRow)]);
                    const int electronRow = phinOffset() + rowNode;
                    const int holeRow = phipOffset() + rowNode;
                    if (derivative != 0.0 &&
                        !constrainedRows[static_cast<std::size_t>(electronRow)]) {
                        triplets.emplace_back(
                            electronRow, columnOffset + columnNode, derivative);
                    }
                    if (derivative != 0.0 &&
                        !constrainedRows[static_cast<std::size_t>(holeRow)]) {
                        triplets.emplace_back(
                            holeRow, columnOffset + columnNode, derivative);
                    }
                }
            }
        }
    }

    SparseMatrixd jacobian(3 * N, 3 * N);
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    return jacobian;
}

SparseMatrixd CoupledDDAssembler::assembleJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs) const
{
    ScopedPerformanceTimer timer("dd.jacobian");
    incrementPerformanceCounter("dd.jacobian_calls");
    const Index Nidx = mesh_.numNodes();
    const int N = static_cast<int>(Nidx);
    if (x.size() != 3 * N)
        throw std::invalid_argument("CoupledDDAssembler::assembleJacobian: vector size mismatch.");

    const VectorXd n = electronDensity(x);
    const VectorXd p = holeDensity(x);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;
    const Real fieldFactor = scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0;
    const VectorXd psi = x.segment(psiOffset(), N) * potentialScale;
    const VectorXd phinState =
        (x.segment(phinOffset(), N) * potentialScale).array()
        + electronQfReference_V_;
    const VectorXd phipState =
        (x.segment(phipOffset(), N) * potentialScale).array()
        + holeQfReference_V_;
    const std::vector<Real> nodeElectricFields =
        (impactIonizationCoupled_ || recombination_.bandToBandEnabled())
        ? detail::computeNodeElectricFields(psi, mesh_, fieldFactor)
        : std::vector<Real>{};
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig_);
    const bool currentAlignedImpact =
        !impactIonizationConfig_.debugRawVanOverstraeten &&
        detail::usesCurrentAlignedAvalancheDrivingForce(impactIonizationConfig_);
    const bool sentaurusEparallelImpact =
        !impactIonizationConfig_.debugRawVanOverstraeten &&
        detail::usesSentaurusEparallelAvalancheDrivingForce(
            impactIonizationConfig_);
    const bool cellCurrentReconstructedCurrent =
        detail::usesCellCurrentReconstructedAvalancheCurrent(impactIonizationConfig_);
    const bool cellVectorCurrentReconstructedCurrent =
        detail::usesCellVectorCurrentReconstructedAvalancheCurrent(impactIonizationConfig_);
    const bool nodalVectorCurrentReconstructedCurrent =
        detail::usesNodalVectorCurrentReconstructedAvalancheCurrent(
            impactIonizationConfig_);
    const bool dualFaceVectorCurrentMagnitude =
        impactIonizationConfig_.currentMagnitudeMode == "dual_face_vector_mag";
    const std::vector<Real> nodeElectronDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phinState, n, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const std::vector<Real> nodeHoleDrivingFields = (impactIonizationCoupled_ && qfImpact)
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig_, mesh_, nodeCells_, psi, phipState, p, ni_, Vt_, fieldFactor)
        : nodeElectricFields;
    const bool qfMobility = mobilityConfig_.highFieldDrivingForce == "quasi_fermi_gradient";
    const bool vectorQfMobility = qfMobility &&
        mobilityConfig_.highFieldGradientDiscretization == "transport_cell_vector";
    const std::vector<Real> electronVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phinState, fieldFactor)
        : std::vector<Real>{};
    const std::vector<Real> holeVectorMobilityFields = vectorQfMobility
        ? detail::transportCellVectorEdgeGradientMagnitudes(
              mesh_, edgeCells_, cellMaterials_, phipState, fieldFactor)
        : std::vector<Real>{};
    const Real chargeAreaFactor = scaling_.enabled ? scaling_.chargeAreaFactor : 1.0;
    const Real sourceIntegralFactor = scaling_.enabled
        ? scaling_.unitSystem.continuitySourceIntegralFactor()
        : 1.0;
    const std::vector<bool> contactNodes = detail::contactNodeMask(mesh_);
    const std::vector<std::vector<Index>> cellEdges = detail::buildCellEdgeMap(edgeCells_, mesh_);
    const std::vector<std::vector<Index>> nodeEdges = detail::buildNodeEdgeMap(mesh_);
    const bool transportMobilityDerivative = transportMobilityDependsOnPotentials(mobilityConfig_);
    const bool sgCurrentAvalanche = impactIonizationCoupled_ &&
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig_);
    const bool triangleGssAvalanche = sgCurrentAvalanche &&
        detail::usesTriangleGssAvalancheSource(impactIonizationConfig_);
    const bool elementEdgeGssLauxAvalanche = sgCurrentAvalanche &&
        detail::usesElementEdgeGssLauxAvalancheSource(
            impactIonizationConfig_);
    const bool cellLocalAvalanche =
        triangleGssAvalanche || elementEdgeGssLauxAvalanche;
    const bool directionalEdgePartition =
        detail::usesDirectionalEdgeAvalancheSourcePartition(impactIonizationConfig_);

    std::vector<bool> constrainedRows(static_cast<std::size_t>(3 * N), false);
    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(psiOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phinOffset() + static_cast<int>(node))] = true;
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        constrainedRows[static_cast<std::size_t>(phipOffset() + static_cast<int>(node))] = true;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(N) * 27);
    std::vector<Real> assembledDiagonal(static_cast<std::size_t>(3 * N), 0.0);
    auto scaleDerivative = [&](int row, int, Real value) {
        if (!scaling_.enabled)
            return value;
        const Real rowScale = (row < N)
            ? (scaling_.permittivityReference_F_per_m * scaling_.V0)
            : (scaling_.C0 * scaling_.D0);
        return value * scaling_.V0 / rowScale;
    };

    auto add = [&](int row, int col, Real value) {
        if (value != 0.0 && !constrainedRows[static_cast<std::size_t>(row)]) {
            const Real scaled = scaleDerivative(row, col, value);
            triplets.emplace_back(row, col, scaled);
            if (row == col)
                assembledDiagonal[static_cast<std::size_t>(row)] += scaled;
        }
    };

    auto addGauge = [&](int row, int col, Real value) {
        if (value != 0.0 && !constrainedRows[static_cast<std::size_t>(row)]) {
            triplets.emplace_back(row, col, value);
            if (row == col)
                assembledDiagonal[static_cast<std::size_t>(row)] += value;
        }
    };

    std::vector<bool> hasElectronContribution(static_cast<std::size_t>(N), false);
    std::vector<bool> hasHoleContribution(static_cast<std::size_t>(N), false);

    struct EdgeAvalancheNodeSources {
        Real node0 = 0.0;
        Real node1 = 0.0;
    };

    // Scharfetter-Gummel edge-current avalanche source for one edge, evaluated
    // directly from the endpoint potentials. This mirrors
    // detail::sgEdgeCurrentAvalancheSourceRecords and returns the configured
    // nodal source partition used by the continuity residuals.
    auto edgeAvalancheNodeSources =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j,
            auto&& phinAt,
            auto&& phipAt) -> EdgeAvalancheNodeSources {
        EdgeAvalancheNodeSources sources;
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return sources;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real niI = ni_[idxI];
        const Real niJ = ni_[idxJ];
        const Real phin_i = phinAt(idxI);
        const Real phin_j = phinAt(idxJ);
        const Real phip_i = phipAt(idxI);
        const Real phip_j = phipAt(idxJ);
        const Real n_i = niI * limitedExp((psi_i - phin_i) / Vt_);
        const Real n_j = niJ * limitedExp((psi_j - phin_j) / Vt_);
        const Real p_i = niI * limitedExp((phip_i - psi_i) / Vt_);
        const Real p_j = niJ * limitedExp((phip_j - psi_j) / Vt_);
        auto psiAt = [&](Index node) {
            const int nodeIndex = static_cast<int>(node);
            return node == idxI ? psi_i : (node == idxJ ? psi_j : psi(nodeIndex));
        };
        auto electronQfAt = [&](Index node) {
            const Real psiNode = psiAt(node);
            const Real phinNode = phinAt(node);
            const Real nNode = ni_[node] * limitedExp((psiNode - phinNode) / Vt_);
            return detail::electronQfForAvalancheGradient(
                psiNode, phinNode, nNode, ni_[node], Vt_, impactIonizationConfig_);
        };
        auto holeQfAt = [&](Index node) {
            const Real psiNode = psiAt(node);
            const Real phipNode = phipAt(node);
            const Real pNode = ni_[node] * limitedExp((phipNode - psiNode) / Vt_);
            return detail::holeQfForAvalancheGradient(
                psiNode, phipNode, pNode, ni_[node], Vt_, impactIonizationConfig_);
        };
        const Real electronQf_i = electronQfAt(idxI);
        const Real electronQf_j = electronQfAt(idxJ);
        const Real holeQf_i = holeQfAt(idxI);
        const Real holeQf_j = holeQfAt(idxJ);
        const Real electricField = std::abs((psi_j - psi_i) / h) * fieldFactor;
        const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h) * fieldFactor;
        const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h) * fieldFactor;
        auto coefficientField = [&](Real edgeQfField, auto&& qfAt) {
            if (detail::usesCellGradientQuasiFermiAvalancheDrive(impactIonizationConfig_)) {
                bool validGradient = false;
                const Point2 gradient = detail::edgeAveragedCellScalarGradient(
                    edgeCells_, mesh_, e, qfAt, validGradient);
                return validGradient ? gradient.norm() : edgeQfField;
            }
            if (impactIonizationConfig_.debugRawVanOverstraeten)
                return edgeQfField;
            return detail::edgeHighFieldDrivingField(
                qfImpact, edgeQfField, electricField, edgeCells_, mesh_, e, contactNodes);
        };
        const Real electronCoefficientField = coefficientField(electronQfField, electronQfAt);
        const Real holeCoefficientField = coefficientField(holeQfField, holeQfAt);
        const Real electronMobilityField = vectorQfMobility
            ? electronVectorMobilityFields[e]
            : (qfMobility ? electronQfField : electricField);
        const Real holeMobilityField = vectorQfMobility
            ? holeVectorMobilityFields[e]
            : (qfMobility ? holeQfField : electricField);
        const Real nAvg = 0.5 * (n_i + n_j);
        const Real pAvg = 0.5 * (p_i + p_j);
        const Real nMid = detail::cellReconstructedAvalancheMidpointDensity(
            impactIonizationConfig_, n_i, n_j, psi_i, psi_j, Vt_);
        const Real pMid = detail::cellReconstructedAvalancheMidpointDensity(
            impactIonizationConfig_, p_i, p_j, psi_j, psi_i, Vt_);
        const Real signedElectricField01 = -(psi_j - psi_i) / h * fieldFactor;
        const Real edgeArea = detail::avalancheSourceEdgeArea(
            impactIonizationConfig_, edgeCells_, mesh_, e, &cellMaterials_);

        auto signedElectronFluxForEdge = [&](Index queryEdge) -> Real {
            const Edge& query = mesh_.getEdge(queryEdge);
            const Real queryH = query.length;
            if (queryH <= 1.0e-30 || couple_[queryEdge] <= 0.0)
                return 0.0;
            const Index a = query.n0;
            const Index b = query.n1;
            const Real psi_a = psiAt(a);
            const Real psi_b = psiAt(b);
            const Real phin_a = phinAt(a);
            const Real phin_b = phinAt(b);
            const Real electronQf_a = electronQfAt(a);
            const Real electronQf_b = electronQfAt(b);
            const Real electric = std::abs((psi_b - psi_a) / queryH);
            const Real mobilityField = qfMobility
                ? std::abs((electronQf_b - electronQf_a) / queryH)
                : electric;
            const Real munLocal = detail::edgeMobility(
                edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, queryEdge,
                CarrierType::Electron, mobilityField, &mobilityConfig_, &psi);
            if (munLocal <= 0.0)
                return 0.0;
            return sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[a], ni_[b], psi_a, psi_b, phin_a, phin_b, Vt_,
                munLocal * Vt_ / queryH, bgnEnabled_);
        };
        auto signedHoleFluxForEdge = [&](Index queryEdge) -> Real {
            const Edge& query = mesh_.getEdge(queryEdge);
            const Real queryH = query.length;
            if (queryH <= 1.0e-30 || couple_[queryEdge] <= 0.0)
                return 0.0;
            const Index a = query.n0;
            const Index b = query.n1;
            const Real psi_a = psiAt(a);
            const Real psi_b = psiAt(b);
            const Real phip_a = phipAt(a);
            const Real phip_b = phipAt(b);
            const Real holeQf_a = holeQfAt(a);
            const Real holeQf_b = holeQfAt(b);
            const Real electric = std::abs((psi_b - psi_a) / queryH);
            const Real mobilityField = qfMobility
                ? std::abs((holeQf_b - holeQf_a) / queryH)
                : electric;
            const Real mupLocal = detail::edgeMobility(
                edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, queryEdge,
                CarrierType::Hole, mobilityField, &mobilityConfig_, &psi);
            if (mupLocal <= 0.0)
                return 0.0;
            return sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[a], ni_[b], psi_a, psi_b, phip_a, phip_b, Vt_,
                mupLocal * Vt_ / queryH, bgnEnabled_);
        };
        auto cellCurrentReconstructedFlux = [&](auto&& signedFluxForEdge) -> Real {
            if (e >= edgeCells_.size())
                return 0.0;
            Real edgeSum = 0.0;
            int adjacentCellCount = 0;
            for (Index cellId : edgeCells_[e]) {
                if (cellId >= cellEdges.size())
                    continue;
                Real cellSum = 0.0;
                int cellEdgeCount = 0;
                for (Index otherEdge : cellEdges[cellId]) {
                    cellSum += std::abs(signedFluxForEdge(otherEdge));
                    ++cellEdgeCount;
                }
                if (cellEdgeCount <= 0)
                    continue;
                edgeSum += cellSum / static_cast<Real>(cellEdgeCount);
                ++adjacentCellCount;
            }
            return adjacentCellCount > 0
                ? edgeSum / static_cast<Real>(adjacentCellCount)
                : 0.0;
        };

        auto cellVectorCurrentReconstructedFlux = [&](auto&& signedFluxForEdge) -> Real {
            if (e >= edgeCells_.size())
                return 0.0;
            Real edgeSum = 0.0;
            int adjacentCellCount = 0;
            for (Index cellId : edgeCells_[e]) {
                if (cellId >= cellEdges.size())
                    continue;
                Real a00 = 0.0;
                Real a01 = 0.0;
                Real a11 = 0.0;
                Real b0 = 0.0;
                Real b1 = 0.0;
                Real absSum = 0.0;
                int used = 0;
                for (Index otherEdge : cellEdges[cellId]) {
                    const Edge& other = mesh_.getEdge(otherEdge);
                    if (other.length <= 1.0e-30)
                        continue;
                    const Node& n0 = mesh_.getNode(other.n0);
                    const Node& n1 = mesh_.getNode(other.n1);
                    const Real tx = (n1.x - n0.x) / other.length;
                    const Real ty = (n1.y - n0.y) / other.length;
                    const Real flux = signedFluxForEdge(otherEdge);
                    a00 += tx * tx;
                    a01 += tx * ty;
                    a11 += ty * ty;
                    b0 += tx * flux;
                    b1 += ty * flux;
                    absSum += std::abs(flux);
                    ++used;
                }
                const Real det = a00 * a11 - a01 * a01;
                const Real scale = std::max({std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0}});
                const Real cellMagnitude = (used >= 2 && std::abs(det) > 1.0e-24 * scale)
                    ? std::sqrt(
                        std::pow((b0 * a11 - b1 * a01) / det, 2) +
                        std::pow((a00 * b1 - a01 * b0) / det, 2))
                    : (used > 0 ? absSum / static_cast<Real>(used) : 0.0);
                if (cellMagnitude <= 0.0)
                    continue;
                edgeSum += cellMagnitude;
                ++adjacentCellCount;
            }
            return adjacentCellCount > 0
                ? edgeSum / static_cast<Real>(adjacentCellCount)
                : 0.0;
        };

        auto dualFaceVectorCurrentReconstructedFlux = [&](auto&& signedFluxForEdge) -> Real {
            std::vector<Real> signedFlux(static_cast<std::size_t>(mesh_.numEdges()), 0.0);
            for (Index otherEdge = 0; otherEdge < mesh_.numEdges(); ++otherEdge)
                signedFlux[static_cast<std::size_t>(otherEdge)] = signedFluxForEdge(otherEdge);
            return detail::medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                e, signedFlux, edgeCells_, cellEdges, mesh_);
        };
        const auto activeTransportEdge = [&](Index queryEdge) {
            if (queryEdge >= edgeCells_.size() ||
                mesh_.getEdge(queryEdge).couple <= 0.0)
                return false;
            for (Index cellId : edgeCells_[queryEdge]) {
                if (cellId < cellMaterials_.size() &&
                    detail::isTransportMaterial(cellMaterials_[cellId]))
                    return true;
            }
            return false;
        };
        auto nodalVectorCurrent = [&](auto&& signedFluxForEdge) -> Point2 {
            return detail::edgeAveragedNodalCurrentVector(
                e, nodeEdges, mesh_, signedFluxForEdge, activeTransportEdge);
        };
        auto nodalVectorCurrentReconstructedFlux =
            [&](auto&& signedFluxForEdge) -> Real {
            return detail::edgeAveragedNodalCurrentMagnitude(
                e, nodeEdges, mesh_, signedFluxForEdge, activeTransportEdge);
        };
        Point2 edgeElectricVector = Point2::Zero();
        if (sentaurusEparallelImpact &&
            nodalVectorCurrentReconstructedCurrent) {
            bool validElectricGradient = false;
            const Point2 electricGradient = detail::edgeAveragedCellScalarGradient(
                edgeCells_, mesh_, e, psiAt, validElectricGradient);
            if (validElectricGradient) {
                edgeElectricVector = -fieldFactor * electricGradient;
            } else {
                const Node& node0 = mesh_.getNode(idxI);
                const Node& node1 = mesh_.getNode(idxJ);
                edgeElectricVector = signedElectricField01 * Point2{
                    (node1.x - node0.x) / h,
                    (node1.y - node0.y) / h};
            }
        }
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, &psi);
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, &psi);
        const Real signedFluxN = mun > 0.0
            ? sgElectronContinuityFluxFromQuasiFermiVariableNi(
                niI, niJ, psi_i, psi_j, phin_i, phin_j, Vt_, mun * Vt_ / h, bgnEnabled_)
            : 0.0;
        const Real signedFluxP = mup > 0.0
            ? sgHoleContinuityFluxFromQuasiFermiVariableNi(
                niI, niJ, psi_i, psi_j, phip_i, phip_j, Vt_, mup * Vt_ / h, bgnEnabled_)
            : 0.0;
        const Real conservedTotalFluxMagnitude =
            detail::conservedTotalCurrentFluxMagnitude(signedFluxN, signedFluxP);

        Real electronSource = 0.0;
        if (mun > 0.0) {
            const Real electronImpactField =
                sentaurusEparallelImpact &&
                nodalVectorCurrentReconstructedCurrent
                ? detail::sentaurusEparallelAvalancheDrivingField(
                    edgeElectricVector,
                    -nodalVectorCurrent(signedElectronFluxForEdge))
                : (currentAlignedImpact
                ? detail::parallelCurrentAvalancheDrivingField(signedElectricField01, signedFluxN)
                : detail::electronAvalancheDrivingField(
                    impactIonizationConfig_, electronCoefficientField, electricField, nAvg));
            const Real rawFluxN = std::abs(signedFluxN);
            const Real reconstructedFluxN = dualFaceVectorCurrentMagnitude
                ? dualFaceVectorCurrentReconstructedFlux(signedElectronFluxForEdge)
                : (nodalVectorCurrentReconstructedCurrent
                    ? nodalVectorCurrentReconstructedFlux(
                        signedElectronFluxForEdge)
                    : (cellVectorCurrentReconstructedCurrent
                    ? cellVectorCurrentReconstructedFlux(signedElectronFluxForEdge)
                    : (cellCurrentReconstructedCurrent
                        ? cellCurrentReconstructedFlux(signedElectronFluxForEdge)
                        : rawFluxN)));
            const Real fluxN = detail::selectAvalancheCurrentFluxProxy(
                impactIonizationConfig_,
                rawFluxN,
                reconstructedFluxN,
                mun,
                nMid,
                electronImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            const Real alphaN = impactIonization_->electronCoefficient(electronImpactField);
            electronSource = alphaN * fluxN * edgeArea;
        }

        Real holeSource = 0.0;
        if (mup > 0.0) {
            const Real holeImpactField =
                sentaurusEparallelImpact &&
                nodalVectorCurrentReconstructedCurrent
                ? detail::sentaurusEparallelAvalancheDrivingField(
                    edgeElectricVector,
                    nodalVectorCurrent(signedHoleFluxForEdge))
                : (currentAlignedImpact
                ? detail::parallelCurrentAvalancheDrivingField(signedElectricField01, signedFluxP)
                : detail::holeAvalancheDrivingField(
                    impactIonizationConfig_, holeCoefficientField, electricField, pAvg));
            const Real rawFluxP = std::abs(signedFluxP);
            const Real reconstructedFluxP = dualFaceVectorCurrentMagnitude
                ? dualFaceVectorCurrentReconstructedFlux(signedHoleFluxForEdge)
                : (nodalVectorCurrentReconstructedCurrent
                    ? nodalVectorCurrentReconstructedFlux(
                        signedHoleFluxForEdge)
                    : (cellVectorCurrentReconstructedCurrent
                    ? cellVectorCurrentReconstructedFlux(signedHoleFluxForEdge)
                    : (cellCurrentReconstructedCurrent
                        ? cellCurrentReconstructedFlux(signedHoleFluxForEdge)
                        : rawFluxP)));
            const Real fluxP = detail::selectAvalancheCurrentFluxProxy(
                impactIonizationConfig_,
                rawFluxP,
                reconstructedFluxP,
                mup,
                pMid,
                holeImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            const Real alphaP = impactIonization_->holeCoefficient(holeImpactField);
            holeSource = alphaP * fluxP * edgeArea;
        }

        detail::EdgeAvalancheDirectionalWeights weights;
        if (directionalEdgePartition) {
            weights = detail::edgeAvalancheDirectionalWeights(
                edgeCells_,
                mesh_,
                e,
                [&](Index node) {
                    const int nodeIndex = static_cast<int>(node);
                    const Real psiNode = node == idxI ? psi_i : (node == idxJ ? psi_j : psi(nodeIndex));
                    const Real phinNode = phinAt(node);
                    const Real nNode = ni_[node] * limitedExp((psiNode - phinNode) / Vt_);
                    return detail::electronQfForAvalancheGradient(
                        psiNode, phinNode, nNode, ni_[node], Vt_, impactIonizationConfig_);
                },
                [&](Index node) {
                    const int nodeIndex = static_cast<int>(node);
                    const Real psiNode = node == idxI ? psi_i : (node == idxJ ? psi_j : psi(nodeIndex));
                    const Real phipNode = phipAt(node);
                    const Real pNode = ni_[node] * limitedExp((phipNode - psiNode) / Vt_);
                    return detail::holeQfForAvalancheGradient(
                        psiNode, phipNode, pNode, ni_[node], Vt_, impactIonizationConfig_);
                });
        }
        sources.node0 = weights.electronNode0 * electronSource + weights.holeNode0 * holeSource;
        sources.node1 = weights.electronNode1 * electronSource + weights.holeNode1 * holeSource;
        return sources;
    };

    auto edgeElectronTransportFlux =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j, Real phin_i, Real phin_j) -> Real {
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return 0.0;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real electronMobilityField = vectorQfMobility
            ? electronVectorMobilityFields[e]
            : (qfMobility
                ? std::abs((phin_j - phin_i) / h) * fieldFactor
                : electricField);
        VectorXd psiForSurface;
        const VectorXd* psiForMobility = &psi;
        if (isSurfaceMobilityModel(mobilityConfig_)) {
            psiForSurface = psi;
            psiForSurface(i) = psi_i;
            psiForSurface(j) = psi_j;
            psiForMobility = &psiForSurface;
        }
        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Electron, electronMobilityField, &mobilityConfig_, psiForMobility);
        if (mun <= 0.0)
            return 0.0;
        const Real coef = mun * Vt_ * fieldFactor * couple_e / h;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real psiRelativeI = psi_i - electronQfReference_V_;
            const Real psiRelativeJ = psi_j - electronQfReference_V_;
            const Real n_i = vela::electronDensity(
                ni_[idxI], Nc_[idxI], psiRelativeI, phin_i, Vt_, carrierStatistics_);
            const Real n_j = vela::electronDensity(
                ni_[idxJ], Nc_[idxJ], psiRelativeJ, phin_j, Vt_, carrierStatistics_);
            const Real etaI = (psiRelativeI - phin_i) / Vt_
                + std::log(ni_[idxI] / Nc_[idxI]);
            const Real etaJ = (psiRelativeJ - phin_j) / Vt_
                + std::log(ni_[idxJ] / Nc_[idxJ]);
            const Real driftPotential = dpsi + Vt_ * std::log(
                (ni_[idxJ] / Nc_[idxJ]) / (ni_[idxI] / Nc_[idxI]));
            return sgElectronFermiDiracContinuityFlux(
                n_i, n_j, etaI, etaJ, driftPotential,
                phin_i, phin_j, Vt_, coef);
        }
        if (bgnEnabled_) {
            return sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ],
                psi_i - electronQfReference_V_,
                psi_j - electronQfReference_V_,
                phin_i, phin_j, Vt_, coef);
        }
        if (ni_[idxI] == ni_[idxJ]) {
            return sgElectronContinuityFluxFromQuasiFermiStable(
                ni_[idxI],
                psi_i - electronQfReference_V_,
                psi_j - electronQfReference_V_,
                phin_i, phin_j, Vt_, coef);
        }
        const Real n_i = ni_[idxI] * limitedExp(
            (psi_i - electronQfReference_V_ - phin_i) / Vt_);
        const Real n_j = ni_[idxJ] * limitedExp(
            (psi_j - electronQfReference_V_ - phin_j) / Vt_);
        return sgElectronContinuityFlux(n_i, n_j, dpsi, Vt_, coef);
    };

    auto edgeHoleTransportFlux =
        [&](Index e, int i, int j, Real h,
            Real psi_i, Real psi_j, Real phip_i, Real phip_j) -> Real {
        const Real couple_e = couple_[e];
        if (h <= 1.0e-30 || couple_e <= 0.0)
            return 0.0;
        const Index idxI = static_cast<Index>(i);
        const Index idxJ = static_cast<Index>(j);
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real holeMobilityField = vectorQfMobility
            ? holeVectorMobilityFields[e]
            : (qfMobility
                ? std::abs((phip_j - phip_i) / h) * fieldFactor
                : electricField);
        VectorXd psiForSurface;
        const VectorXd* psiForMobility = &psi;
        if (isSurfaceMobilityModel(mobilityConfig_)) {
            psiForSurface = psi;
            psiForSurface(i) = psi_i;
            psiForSurface(j) = psi_j;
            psiForMobility = &psiForSurface;
        }
        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e,
            CarrierType::Hole, holeMobilityField, &mobilityConfig_, psiForMobility);
        if (mup <= 0.0)
            return 0.0;
        const Real coef = mup * Vt_ * fieldFactor * couple_e / h;
        if (usesFermiDirac(carrierStatistics_)) {
            const Real psiRelativeI = psi_i - holeQfReference_V_;
            const Real psiRelativeJ = psi_j - holeQfReference_V_;
            const Real p_i = vela::holeDensity(
                ni_[idxI], Nv_[idxI], psiRelativeI, phip_i, Vt_, carrierStatistics_);
            const Real p_j = vela::holeDensity(
                ni_[idxJ], Nv_[idxJ], psiRelativeJ, phip_j, Vt_, carrierStatistics_);
            const Real etaI = (phip_i - psiRelativeI) / Vt_
                + std::log(ni_[idxI] / Nv_[idxI]);
            const Real etaJ = (phip_j - psiRelativeJ) / Vt_
                + std::log(ni_[idxJ] / Nv_[idxJ]);
            const Real driftPotential = dpsi + Vt_ * std::log(
                (ni_[idxI] / Nv_[idxI]) / (ni_[idxJ] / Nv_[idxJ]));
            return sgHoleFermiDiracContinuityFlux(
                p_i, p_j, etaI, etaJ, driftPotential,
                phip_i, phip_j, Vt_, coef);
        }
        if (bgnEnabled_) {
            return sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni_[idxI], ni_[idxJ],
                psi_i - holeQfReference_V_,
                psi_j - holeQfReference_V_,
                phip_i, phip_j, Vt_, coef);
        }
        if (ni_[idxI] == ni_[idxJ]) {
            return sgHoleContinuityFluxFromQuasiFermiStable(
                ni_[idxI],
                psi_i - holeQfReference_V_,
                psi_j - holeQfReference_V_,
                phip_i, phip_j, Vt_, coef);
        }
        const Real p_i = ni_[idxI] * limitedExp(
            (phip_i - (psi_i - holeQfReference_V_)) / Vt_);
        const Real p_j = ni_[idxJ] * limitedExp(
            (phip_j - (psi_j - holeQfReference_V_)) / Vt_);
        return sgHoleContinuityFlux(p_i, p_j, dpsi, Vt_, coef);
    };

    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real h = edge.length;
        if (h < 1.0e-30) continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = x(psiOffset() + i) * potentialScale;
        const Real psi_j = x(psiOffset() + j) * potentialScale;
        const Real phin_i = x(phinOffset() + i) * potentialScale;
        const Real phin_j = x(phinOffset() + j) * potentialScale;
        const Real phip_i = x(phipOffset() + i) * potentialScale;
        const Real phip_j = x(phipOffset() + j) * potentialScale;
        const Real dpsi = psi_j - psi_i;
        const Real electricField = std::abs(dpsi / h) * fieldFactor;
        const Real electronMobilityField =
            vectorQfMobility ? electronVectorMobilityFields[e] :
            (qfMobility ? std::abs((phin_j - phin_i) / h) * fieldFactor : electricField);
        const Real holeMobilityField =
            vectorQfMobility ? holeVectorMobilityFields[e] :
            (qfMobility ? std::abs((phip_j - phip_i) / h) * fieldFactor : electricField);
        const Real u = dpsi / Vt_;
        const Real Bu = bernoulli(u);
        const Real dBu = bernoulliDerivative(u);
        const Real Bminus = bernoulli(-u);
        const Real dBminusDu = -bernoulliDerivative(-u);

        const Real eps = detail::edgeEpsilon(edgeCells_, mesh_, matdb_, e);
        const Real G = eps * couple_[e] / h;
        add(psiOffset() + i, psiOffset() + i,  G);
        add(psiOffset() + i, psiOffset() + j, -G);
        add(psiOffset() + j, psiOffset() + i, -G);
        add(psiOffset() + j, psiOffset() + j,  G);

        const Real mun = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Electron,
            electronMobilityField,
            &mobilityConfig_,
            &psi);
        if (mun > 0.0) {
            hasElectronContribution[static_cast<std::size_t>(i)] = true;
            hasElectronContribution[static_cast<std::size_t>(j)] = true;

            if (transportMobilityDerivative || usesFermiDirac(carrierStatistics_)) {
                const Real vals[4] = {psi_i, psi_j, phin_i, phin_j};
                const int cols[4] = {
                    psiOffset() + i, psiOffset() + j,
                    phinOffset() + i, phinOffset() + j,
                };
                for (int k = 0; k < 4; ++k) {
                    const Real step = 1.0e-6 * std::max(1.0, std::abs(vals[k]));
                    Real vp[4] = {vals[0], vals[1], vals[2], vals[3]};
                    Real vm[4] = {vals[0], vals[1], vals[2], vals[3]};
                    vp[k] += step;
                    vm[k] -= step;
                    const Real fp = edgeElectronTransportFlux(
                        e, i, j, h, vp[0], vp[1], vp[2], vp[3]);
                    const Real fm = edgeElectronTransportFlux(
                        e, i, j, h, vm[0], vm[1], vm[2], vm[3]);
                    const Real dF = (fp - fm) / (2.0 * step);
                    add(phinOffset() + i, cols[k], dF);
                    add(phinOffset() + j, cols[k], -dF);
                }
            } else {
                const Real coef = mun * Vt_ * fieldFactor * couple_[e] / h;
                Real dF_dpsi_i = 0.0;
                Real dF_dpsi_j = 0.0;
                Real dF_dphin_i = 0.0;
                Real dF_dphin_j = 0.0;
                const Index idxI = static_cast<Index>(i);
                const Index idxJ = static_cast<Index>(j);
                const Real eta = u + std::log(ni_[idxJ] / ni_[idxI]);
                const Real Bplus = bernoulli(eta);
                const Real Bminus = bernoulli(-eta);
                const Real dBplus = bernoulliDerivative(eta);
                const Real dBminusArg = bernoulliDerivative(-eta);
                const Real niI = ni_[idxI];
                const Real niJ = ni_[idxJ];
                const Real nI = niI * limitedExp(
                    (psi_i - electronQfReference_V_ - phin_i) / Vt_);
                const Real nJ = niJ * limitedExp(
                    (psi_j - electronQfReference_V_ - phin_j) / Vt_);
                dF_dpsi_i = coef / Vt_ * ((dBminusArg + Bminus) * nI + dBplus * nJ);
                dF_dpsi_j = coef / Vt_ * (-dBminusArg * nI - (dBplus + Bplus) * nJ);
                dF_dphin_i = coef * Bminus * (-nI / Vt_);
                dF_dphin_j = coef * Bplus * ( nJ / Vt_);

                add(phinOffset() + i, psiOffset() + i, dF_dpsi_i);
                add(phinOffset() + i, psiOffset() + j, dF_dpsi_j);
                add(phinOffset() + i, phinOffset() + i, dF_dphin_i);
                add(phinOffset() + i, phinOffset() + j, dF_dphin_j);
                add(phinOffset() + j, psiOffset() + i, -dF_dpsi_i);
                add(phinOffset() + j, psiOffset() + j, -dF_dpsi_j);
                add(phinOffset() + j, phinOffset() + i, -dF_dphin_i);
                add(phinOffset() + j, phinOffset() + j, -dF_dphin_j);
            }
        }

        const Real mup = detail::edgeMobility(
            edgeCells_, mesh_, doping_, *mobility_, cellMaterials_, e, CarrierType::Hole,
            holeMobilityField,
            &mobilityConfig_,
            &psi);
        if (mup > 0.0) {
            hasHoleContribution[static_cast<std::size_t>(i)] = true;
            hasHoleContribution[static_cast<std::size_t>(j)] = true;

            if (transportMobilityDerivative || usesFermiDirac(carrierStatistics_)) {
                const Real vals[4] = {psi_i, psi_j, phip_i, phip_j};
                const int cols[4] = {
                    psiOffset() + i, psiOffset() + j,
                    phipOffset() + i, phipOffset() + j,
                };
                for (int k = 0; k < 4; ++k) {
                    const Real step = 1.0e-6 * std::max(1.0, std::abs(vals[k]));
                    Real vp[4] = {vals[0], vals[1], vals[2], vals[3]};
                    Real vm[4] = {vals[0], vals[1], vals[2], vals[3]};
                    vp[k] += step;
                    vm[k] -= step;
                    const Real fp = edgeHoleTransportFlux(
                        e, i, j, h, vp[0], vp[1], vp[2], vp[3]);
                    const Real fm = edgeHoleTransportFlux(
                        e, i, j, h, vm[0], vm[1], vm[2], vm[3]);
                    const Real dF = (fp - fm) / (2.0 * step);
                    add(phipOffset() + i, cols[k], dF);
                    add(phipOffset() + j, cols[k], -dF);
                }
            } else {
                const Real coef = mup * Vt_ * fieldFactor * couple_[e] / h;
                Real dF_dpsi_i = 0.0;
                Real dF_dpsi_j = 0.0;
                Real dF_dphip_i = 0.0;
                Real dF_dphip_j = 0.0;
                const Index idxI = static_cast<Index>(i);
                const Index idxJ = static_cast<Index>(j);
                const Real eta = u + std::log(ni_[idxI] / ni_[idxJ]);
                const Real Bplus = bernoulli(eta);
                const Real Bminus = bernoulli(-eta);
                const Real dBplus = bernoulliDerivative(eta);
                const Real dBminusArg = bernoulliDerivative(-eta);
                const Real niI = ni_[idxI];
                const Real niJ = ni_[idxJ];
                const Real pI = niI * limitedExp(
                    (phip_i - (psi_i - holeQfReference_V_)) / Vt_);
                const Real pJ = niJ * limitedExp(
                    (phip_j - (psi_j - holeQfReference_V_)) / Vt_);
                dF_dpsi_i = coef / Vt_ * (-(dBplus + Bplus) * pI - dBminusArg * pJ);
                dF_dpsi_j = coef / Vt_ * (dBplus * pI + (dBminusArg + Bminus) * pJ);
                dF_dphip_i = coef * Bplus * ( pI / Vt_);
                dF_dphip_j = coef * Bminus * (-pJ / Vt_);

                add(phipOffset() + i, psiOffset() + i, dF_dpsi_i);
                add(phipOffset() + i, psiOffset() + j, dF_dpsi_j);
                add(phipOffset() + i, phipOffset() + i, dF_dphip_i);
                add(phipOffset() + i, phipOffset() + j, dF_dphip_j);
                add(phipOffset() + j, psiOffset() + i, -dF_dpsi_i);
                add(phipOffset() + j, psiOffset() + j, -dF_dpsi_j);
                add(phipOffset() + j, phipOffset() + i, -dF_dphip_i);
                add(phipOffset() + j, phipOffset() + j, -dF_dphip_j);
            }
        }

        if (sgCurrentAvalanche && !cellLocalAvalanche) {
            // Finite-difference the directionally partitioned nodal avalanche
            // source with respect to the six endpoint potentials. This captures
            // the flux (carrier density), driving-field (alpha), field-dependent
            // mobility, and qF-gradient partition derivatives together.
            auto phinAt = [&](Index node) { return phinState(static_cast<int>(node)); };
            auto phipAt = [&](Index node) { return phipState(static_cast<int>(node)); };
            const EdgeAvalancheNodeSources base = edgeAvalancheNodeSources(
                e, i, j, h, psi_i, psi_j, phinAt, phipAt);
            std::vector<Index> qfStencilNodes = {static_cast<Index>(i), static_cast<Index>(j)};
            if (e < edgeCells_.size()) {
                for (const Index cellId : edgeCells_[e]) {
                    const Cell& cell = mesh_.getCell(cellId);
                    for (const Index node : cell.node_ids) {
                        if (std::find(qfStencilNodes.begin(), qfStencilNodes.end(), node) ==
                            qfStencilNodes.end()) {
                            qfStencilNodes.push_back(node);
                        }
                    }
                }
            }

            std::vector<int> cols;
            std::vector<Real> dS0;
            std::vector<Real> dS1;
            bool anyNonzero = false;
            auto appendDerivative = [&](int col, const EdgeAvalancheNodeSources& sp,
                                        const EdgeAvalancheNodeSources& sm, Real step) {
                cols.push_back(col);
                dS0.push_back(
                    (sp.node0 - sm.node0) * sourceIntegralFactor / (2.0 * step));
                dS1.push_back(
                    (sp.node1 - sm.node1) * sourceIntegralFactor / (2.0 * step));
                if (dS0.back() != 0.0 || dS1.back() != 0.0)
                    anyNonzero = true;
            };

            const Real psiVals[2] = {psi_i, psi_j};
            const int psiCols[2] = {psiOffset() + i, psiOffset() + j};
            for (int k = 0; k < 2; ++k) {
                const Real step = 1.0e-7 * std::max(1.0, std::abs(psiVals[k]));
                Real psiP[2] = {psi_i, psi_j};
                Real psiM[2] = {psi_i, psi_j};
                psiP[k] += step;
                psiM[k] -= step;
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psiP[0], psiP[1], phinAt, phipAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psiM[0], psiM[1], phinAt, phipAt);
                appendDerivative(psiCols[k], sp, sm, step);
            }

            for (const Index node : qfStencilNodes) {
                const int nodeIndex = static_cast<int>(node);
                const Real phinValue = phinState(nodeIndex);
                const Real step = 1.0e-7 * std::max(1.0, std::abs(phinValue));
                auto phinPlusAt = [&](Index queryNode) {
                    return queryNode == node ? phinValue + step
                                             : phinState(static_cast<int>(queryNode));
                };
                auto phinMinusAt = [&](Index queryNode) {
                    return queryNode == node ? phinValue - step
                                             : phinState(static_cast<int>(queryNode));
                };
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinPlusAt, phipAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinMinusAt, phipAt);
                appendDerivative(phinOffset() + nodeIndex, sp, sm, step);
            }

            for (const Index node : qfStencilNodes) {
                const int nodeIndex = static_cast<int>(node);
                const Real phipValue = phipState(nodeIndex);
                const Real step = 1.0e-7 * std::max(1.0, std::abs(phipValue));
                auto phipPlusAt = [&](Index queryNode) {
                    return queryNode == node ? phipValue + step
                                             : phipState(static_cast<int>(queryNode));
                };
                auto phipMinusAt = [&](Index queryNode) {
                    return queryNode == node ? phipValue - step
                                             : phipState(static_cast<int>(queryNode));
                };
                const EdgeAvalancheNodeSources sp = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinAt, phipPlusAt);
                const EdgeAvalancheNodeSources sm = edgeAvalancheNodeSources(
                    e, i, j, h, psi_i, psi_j, phinAt, phipMinusAt);
                appendDerivative(phipOffset() + nodeIndex, sp, sm, step);
            }

            if (base.node0 != 0.0 || base.node1 != 0.0 || anyNonzero) {
                const int node0Rows[2] = {phinOffset() + i, phipOffset() + i};
                const int node1Rows[2] = {phinOffset() + j, phipOffset() + j};
                for (int row : node0Rows)
                    for (std::size_t k = 0; k < cols.size(); ++k)
                        add(row, cols[k], -dS0[k]);
                for (int row : node1Rows)
                    for (std::size_t k = 0; k < cols.size(); ++k)
                        add(row, cols[k], -dS1[k]);
                hasElectronContribution[static_cast<std::size_t>(i)] = true;
                hasElectronContribution[static_cast<std::size_t>(j)] = true;
                hasHoleContribution[static_cast<std::size_t>(i)] = true;
                hasHoleContribution[static_cast<std::size_t>(j)] = true;
            }
        }
    }

    if (cellLocalAvalanche) {
        if (isSurfaceMobilityModel(mobilityConfig_)) {
            throw std::invalid_argument(
                "cell-local avalanche source does not support surface mobility");
        }
        if (elementEdgeGssLauxAvalanche) {
            for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
                const Cell& cell = mesh_.getCell(cellId);
                std::array<detail::Tri3LocalForwardDual, 3> localPsi{};
                std::array<detail::Tri3LocalForwardDual, 3> localPhin{};
                std::array<detail::Tri3LocalForwardDual, 3> localPhip{};
                std::array<detail::Tri3LocalForwardDual, 3> localElectronDensity{};
                std::array<detail::Tri3LocalForwardDual, 3> localHoleDensity{};
                std::array<Real, 3> localIntrinsicDensity{};
                for (std::size_t localNode = 0; localNode < 3; ++localNode) {
                    const Index node = cell.node_ids[localNode];
                    const int nodeIndex = static_cast<int>(node);
                    localPsi[localNode] = detail::Tri3LocalForwardDual::variable(
                        psi(nodeIndex), localNode);
                    localPhin[localNode] = detail::Tri3LocalForwardDual::variable(
                        phinState(nodeIndex), 3 + localNode);
                    localPhip[localNode] = detail::Tri3LocalForwardDual::variable(
                        phipState(nodeIndex), 6 + localNode);
                    localIntrinsicDensity[localNode] = ni_[node];
                    localElectronDensity[localNode] =
                        detail::Tri3LocalForwardDual(ni_[node]) *
                        detail::localAdLimitedExp(
                            (localPsi[localNode] - localPhin[localNode]) /
                            detail::Tri3LocalForwardDual(Vt_));
                    localHoleDensity[localNode] =
                        detail::Tri3LocalForwardDual(ni_[node]) *
                        detail::localAdLimitedExp(
                            (localPhip[localNode] - localPsi[localNode]) /
                            detail::Tri3LocalForwardDual(Vt_));
                }
                const auto sources =
                    detail::elementEdgeGssLauxAvalancheSourceIntegralsLocal<
                        detail::Tri3LocalForwardDual>(
                        impactIonizationConfig_, mobilityConfig_, *mobility_,
                        cellEdges.at(static_cast<std::size_t>(cellId)),
                        mesh_, doping_, cellMaterials_, cellId, localPsi,
                        localPhin, localPhip, localElectronDensity,
                        localHoleDensity, localIntrinsicDensity, Vt_, fieldFactor);

                bool anySourceOrDerivative = false;
                for (int localRow = 0; localRow < 3; ++localRow) {
                    const auto& source =
                        sources.combined[static_cast<std::size_t>(localRow)];
                    anySourceOrDerivative = anySourceOrDerivative || source.value != 0.0;
                    const int rowNode = static_cast<int>(
                        cell.node_ids[static_cast<std::size_t>(localRow)]);
                    for (std::size_t localDof = 0;
                         localDof < detail::Tri3LocalPotentialDofCount; ++localDof) {
                        const int variableBlock = static_cast<int>(localDof / 3);
                        const int localColumn = static_cast<int>(localDof % 3);
                        const int columnNode = static_cast<int>(
                            cell.node_ids[static_cast<std::size_t>(localColumn)]);
                        const int columnOffset = variableBlock == 0
                            ? psiOffset()
                            : (variableBlock == 1 ? phinOffset() : phipOffset());
                        const Real derivative =
                            source.derivative[localDof] * sourceIntegralFactor;
                        anySourceOrDerivative =
                            anySourceOrDerivative || derivative != 0.0;
                        add(phinOffset() + rowNode, columnOffset + columnNode,
                            -derivative);
                        add(phipOffset() + rowNode, columnOffset + columnNode,
                            -derivative);
                    }
                }
                if (anySourceOrDerivative) {
                    for (const Index node : cell.node_ids) {
                        hasElectronContribution[static_cast<std::size_t>(node)] = true;
                        hasHoleContribution[static_cast<std::size_t>(node)] = true;
                    }
                }
            }
        } else {
            auto cellNodeSources = [&](Index cellId,
                                       const VectorXd& psiValues,
                                       const VectorXd& phinValues,
                                       const VectorXd& phipValues,
                                       const VectorXd& nValues,
                                       const VectorXd& pValues) {
                VectorXd sources = VectorXd::Zero(3);
                const Cell& cell = mesh_.getCell(cellId);
                if (elementEdgeGssLauxAvalanche) {
                    const auto record =
                        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
                            impactIonizationConfig_, *impactIonization_,
                            mobilityConfig_, *mobility_,
                            cellEdges.at(static_cast<std::size_t>(cellId)),
                            mesh_, doping_, cellMaterials_, cellId, psiValues,
                            phinValues, phipValues, nValues, pValues, ni_, Vt_,
                            fieldFactor);
                    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
                        sources(static_cast<int>(localNode)) +=
                            record.combinedSourceIntegrals[localNode];
                    }
                    return sources;
                }
                const auto records = detail::triangleGssAvalancheSourceRecordsForCell(
                    impactIonizationConfig_, *impactIonization_, mobilityConfig_, *mobility_,
                    cellEdges.at(static_cast<std::size_t>(cellId)),
                    mesh_, doping_, cellMaterials_, cellId, psiValues, phinValues,
                    phipValues, nValues, pValues, Vt_, fieldFactor);
                for (const auto& record : records) {
                    for (int localNode = 0; localNode < 3; ++localNode) {
                        const Index node =
                            cell.node_ids[static_cast<std::size_t>(localNode)];
                        if (node == record.node0)
                            sources(localNode) += record.node0SourceIntegral;
                        if (node == record.node1)
                            sources(localNode) += record.node1SourceIntegral;
                    }
                }
                return sources;
            };

            for (Index cellId = 0; cellId < mesh_.numCells(); ++cellId) {
                const Cell& cell = mesh_.getCell(cellId);
                const VectorXd baseSources =
                    cellNodeSources(cellId, psi, phinState, phipState, n, p);

                for (int variableBlock = 0; variableBlock < 3; ++variableBlock) {
                    for (int localColumn = 0; localColumn < 3; ++localColumn) {
                        const Index node =
                            cell.node_ids[static_cast<std::size_t>(localColumn)];
                        const int nodeIndex = static_cast<int>(node);
                        const Real value = variableBlock == 0 ? psi(nodeIndex)
                            : (variableBlock == 1 ? phinState(nodeIndex)
                                                  : phipState(nodeIndex));
                        // Match the assembler-level relative perturbation in scaled
                        // coordinates: dV = relativeStep * max(V0, |V|).
                        const Real step =
                            detail::physicalPotentialCentralDifferenceStep(
                                value, potentialScale, 1.0e-7);
                        VectorXd psiPlus = psi;
                        VectorXd psiMinus = psi;
                        VectorXd phinPlus = phinState;
                        VectorXd phinMinus = phinState;
                        VectorXd phipPlus = phipState;
                        VectorXd phipMinus = phipState;
                        if (variableBlock == 0) {
                            psiPlus(nodeIndex) += step;
                            psiMinus(nodeIndex) -= step;
                        } else if (variableBlock == 1) {
                            phinPlus(nodeIndex) += step;
                            phinMinus(nodeIndex) -= step;
                        } else {
                            phipPlus(nodeIndex) += step;
                            phipMinus(nodeIndex) -= step;
                        }

                        VectorXd nPlus = n;
                        VectorXd nMinus = n;
                        VectorXd pPlus = p;
                        VectorXd pMinus = p;
                        nPlus(nodeIndex) = ni_[node] * limitedExp(
                            (psiPlus(nodeIndex) - phinPlus(nodeIndex)) / Vt_);
                        nMinus(nodeIndex) = ni_[node] * limitedExp(
                            (psiMinus(nodeIndex) - phinMinus(nodeIndex)) / Vt_);
                        pPlus(nodeIndex) = ni_[node] * limitedExp(
                            (phipPlus(nodeIndex) - psiPlus(nodeIndex)) / Vt_);
                        pMinus(nodeIndex) = ni_[node] * limitedExp(
                            (phipMinus(nodeIndex) - psiMinus(nodeIndex)) / Vt_);
                        const VectorXd plusSources = cellNodeSources(
                            cellId, psiPlus, phinPlus, phipPlus, nPlus, pPlus);
                        const VectorXd minusSources = cellNodeSources(
                            cellId, psiMinus, phinMinus, phipMinus, nMinus, pMinus);
                        const VectorXd derivative =
                            (plusSources - minusSources) * sourceIntegralFactor /
                            (2.0 * step);
                        const int column = (variableBlock == 0 ? psiOffset()
                            : (variableBlock == 1 ? phinOffset() : phipOffset())) +
                            nodeIndex;
                        for (int localRow = 0; localRow < 3; ++localRow) {
                            const int rowNode = static_cast<int>(
                                cell.node_ids[static_cast<std::size_t>(localRow)]);
                            add(phinOffset() + rowNode, column, -derivative(localRow));
                            add(phipOffset() + rowNode, column, -derivative(localRow));
                        }
                    }
                }
                if (baseSources.squaredNorm() > 0.0) {
                    for (const Index node : cell.node_ids) {
                        hasElectronContribution[static_cast<std::size_t>(node)] = true;
                        hasHoleContribution[static_cast<std::size_t>(node)] = true;
                    }
                }
            }
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        const Real psiRelativeN = psi(ii) - electronQfReference_V_;
        const Real psiRelativeP = psi(ii) - holeQfReference_V_;
        const Real dnDeta = electronDensityDerivativeEta(
            ni_[i], Nc_[i], psiRelativeN, phinState(ii), Vt_, carrierStatistics_);
        const Real dpDeta = holeDensityDerivativeEta(
            ni_[i], Nv_[i], psiRelativeP, phipState(ii), Vt_, carrierStatistics_);
        const Real dni_dpsi = dnDeta / Vt_;
        const Real dni_dphin = -dnDeta / Vt_;
        const Real dpi_dpsi = -dpDeta / Vt_;
        const Real dpi_dphip = dpDeta / Vt_;

        add(psiOffset() + ii, psiOffset() + ii,
            -constants::q * (dpi_dpsi - dni_dpsi) * vol_[i] * chargeAreaFactor);
        add(psiOffset() + ii, phinOffset() + ii,
            -constants::q * (-dni_dphin) * vol_[i] * chargeAreaFactor);
        add(psiOffset() + ii, phipOffset() + ii,
            -constants::q * dpi_dphip * vol_[i] * chargeAreaFactor);

        const Real ni = ni_[i];
        if (ni <= 0.0)
            continue;

        if (recombination_.srhEnabled() || recombination_.augerEnabled()) {
            if (usesFermiDirac(carrierStatistics_)) {
                auto localRate = [&](Real psiValue, Real phinValue, Real phipValue) {
                    const Real nValue = vela::electronDensity(
                        ni, Nc_[i], psiValue - electronQfReference_V_,
                        phinValue - electronQfReference_V_, Vt_, carrierStatistics_);
                    const Real pValue = vela::holeDensity(
                        ni, Nv_[i], psiValue - holeQfReference_V_,
                        phipValue - holeQfReference_V_, Vt_, carrierStatistics_);
                    if (phinValue == phipValue)
                        return Real{0.0};
                    const Real equilibriumProduct = equilibriumCarrierProduct(
                        nValue, pValue, ni, Nc_[i], Nv_[i], Vt_, carrierStatistics_);
                    return recombination_.totalRateFromExcessProduct(
                        nValue * pValue - equilibriumProduct,
                        nValue, pValue,
                        std::sqrt(std::max<Real>(equilibriumProduct, 0.0)));
                };
                const Real values[3] = {psi(ii), phinState(ii), phipState(ii)};
                const int columns[3] = {
                    psiOffset() + ii, phinOffset() + ii, phipOffset() + ii};
                for (int k = 0; k < 3; ++k) {
                    const Real step = 1.0e-7 * std::max<Real>(1.0, std::abs(values[k]));
                    Real plus[3] = {values[0], values[1], values[2]};
                    Real minus[3] = {values[0], values[1], values[2]};
                    plus[k] += step;
                    minus[k] -= step;
                    const Real derivative = (
                        localRate(plus[0], plus[1], plus[2])
                        - localRate(minus[0], minus[1], minus[2])) / (2.0 * step);
                    add(phinOffset() + ii, columns[k],
                        derivative * vol_[i] * sourceIntegralFactor);
                    add(phipOffset() + ii, columns[k],
                        derivative * vol_[i] * sourceIntegralFactor);
                }
                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            } else {
                const Real dPhi =
                    (holeQfReference_V_ - electronQfReference_V_)
                    + (x(phipOffset() + ii) - x(phinOffset() + ii)) * potentialScale;
                const Real np = ni * ni * std::exp(dPhi / Vt_);
                const Real excessProduct = ni * ni * std::expm1(dPhi / Vt_);
                const Real R = recombination_.totalRateFromExcessProduct(
                    excessProduct, n(ii), p(ii), ni);
                if (R != 0.0) {
                const auto deriv = recombination_.totalRateDerivativesFromExcessProduct(
                    excessProduct, n(ii), p(ii), ni);

                const Real dExcess_dphin = -np / Vt_;
                const Real dExcess_dphip =  np / Vt_;
                const Real dR_dpsi = deriv.dRateDn * dni_dpsi + deriv.dRateDp * dpi_dpsi;
                const Real dR_dphin = deriv.dRateDn * dni_dphin
                                     + deriv.dRateDExcess * dExcess_dphin;
                const Real dR_dphip = deriv.dRateDp * dpi_dphip
                                     + deriv.dRateDExcess * dExcess_dphip;

                add(phinOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i] * sourceIntegralFactor);
                add(phinOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i] * sourceIntegralFactor);
                add(phinOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, psiOffset() + ii, dR_dpsi * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, phinOffset() + ii, dR_dphin * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, phipOffset() + ii, dR_dphip * vol_[i] * sourceIntegralFactor);

                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
                }
            }
        }

        if (impactIonizationCoupled_ && sgCurrentAvalanche) {
            // The SG edge-current avalanche probe has nonlocal edge derivatives.
            // Omit source derivatives here rather than applying the legacy
            // node-local approximation to the wrong discretization.
        } else if (impactIonizationCoupled_) {
            const Real electronImpactField = detail::electronAvalancheDrivingField(
                impactIonizationConfig_, nodeElectronDrivingFields[i], nodeElectricFields[i], n(ii));
            const Real holeImpactField = detail::holeAvalancheDrivingField(
                impactIonizationConfig_, nodeHoleDrivingFields[i], nodeElectricFields[i], p(ii));
            const Real alphaN = impactIonization_->electronCoefficient(electronImpactField);
            const Real alphaP = impactIonization_->holeCoefficient(holeImpactField);
            const Real G = detail::impactIonizationGenerationRate(
                impactIonizationConfig_,
                *impactIonization_,
                mobilityConfig_,
                *mobility_,
                nodeCells_,
                mesh_,
                doping_,
                cellMaterials_,
                i,
                nodeElectricFields[i],
                nodeElectronDrivingFields[i],
                nodeHoleDrivingFields[i],
                n(ii),
                p(ii));
            if (G != 0.0) {
                // Local carrier-density derivatives are included in the analytic Jacobian.
                // Driving-field and mobility derivatives are intentionally omitted because
                // both use edge/node max operations; finite-difference Jacobian remains
                // available for exact derivatives of configured avalanche runs.
                Real electronFactor = 0.0;
                Real holeFactor = 0.0;
                if (impactIonizationConfig_.generation == "current_density") {
                    const Real mun = detail::nodeMobility(
                        nodeCells_,
                        mesh_,
                        doping_,
                        *mobility_,
                        cellMaterials_,
                        i,
                        CarrierType::Electron,
                        electronImpactField,
                        &mobilityConfig_);
                    const Real mup = detail::nodeMobility(
                        nodeCells_,
                        mesh_,
                        doping_,
                        *mobility_,
                        cellMaterials_,
                        i,
                        CarrierType::Hole,
                        holeImpactField,
                        &mobilityConfig_);
                    electronFactor = alphaN * mun * std::abs(electronImpactField);
                    holeFactor = alphaP * mup * std::abs(holeImpactField);
                } else {
                    const Real denominator = alphaN * n(ii) + alphaP * p(ii);
                    const Real velocity = (denominator != 0.0) ? (G / denominator) : 0.0;
                    electronFactor = velocity * alphaN;
                    holeFactor = velocity * alphaP;
                }
                const Real dG_dpsi = electronFactor * dni_dpsi + holeFactor * dpi_dpsi;
                const Real dG_dphin = electronFactor * dni_dphin;
                const Real dG_dphip = holeFactor * dpi_dphip;

                add(phinOffset() + ii, psiOffset() + ii, -dG_dpsi * vol_[i] * sourceIntegralFactor);
                add(phinOffset() + ii, phinOffset() + ii, -dG_dphin * vol_[i] * sourceIntegralFactor);
                add(phinOffset() + ii, phipOffset() + ii, -dG_dphip * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, psiOffset() + ii, -dG_dpsi * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, phinOffset() + ii, -dG_dphin * vol_[i] * sourceIntegralFactor);
                add(phipOffset() + ii, phipOffset() + ii, -dG_dphip * vol_[i] * sourceIntegralFactor);

                hasElectronContribution[static_cast<std::size_t>(ii)] = true;
                hasHoleContribution[static_cast<std::size_t>(ii)] = true;
            }
        }
    }

    for (Index i = 0; i < Nidx; ++i) {
        const int ii = static_cast<int>(i);
        if (!hasElectronContribution[static_cast<std::size_t>(ii)])
            addGauge(phinOffset() + ii, phinOffset() + ii, 1.0);
        if (!hasHoleContribution[static_cast<std::size_t>(ii)])
            addGauge(phipOffset() + ii, phipOffset() + ii, 1.0);
    }

    if (recombination_.bandToBandEnabled() &&
        recombination_.bandToBand().config().jacobian ==
            "potential_finite_difference") {
        const Real relativeStep =
            recombination_.bandToBand().config().jacobianRelativeStep;
        auto sourcesForPotential = [&](const VectorXd& potential) {
            return detail::bandToBandGenerationNodeSourceIntegrals(
                recombination_.bandToBand(), scaling_.unitSystem, mesh_,
                cellMaterials_, potential, fieldFactor);
        };

        VectorXd plusPotential = psi;
        VectorXd minusPotential = psi;
        for (Index columnNode = 0; columnNode < Nidx; ++columnNode) {
            const int column = static_cast<int>(columnNode);
            const Real step = relativeStep *
                std::max<Real>(1.0, std::abs(psi(column)));
            plusPotential(column) += step;
            minusPotential(column) -= step;
            const std::vector<Real> plusSources =
                sourcesForPotential(plusPotential);
            const std::vector<Real> minusSources =
                sourcesForPotential(minusPotential);
            plusPotential(column) = psi(column);
            minusPotential(column) = psi(column);

            for (Index rowNode = 0; rowNode < Nidx; ++rowNode) {
                if (ni_[rowNode] <= 0.0)
                    continue;
                const Real derivative =
                    (plusSources[rowNode] - minusSources[rowNode]) / (2.0 * step);
                if (derivative == 0.0)
                    continue;
                const int row = static_cast<int>(rowNode);
                const Real contribution = -derivative * sourceIntegralFactor;
                add(phinOffset() + row, psiOffset() + column, contribution);
                add(phipOffset() + row, psiOffset() + column, contribution);
            }
        }
    }

    if (carrierDiagonalFloor_.enabled &&
        carrierDiagonalFloor_.scale > 0.0 &&
        recombination_.srhEnabled()) {
        const Real tauSum = recombinationConfig_.taun + recombinationConfig_.taup;
        if (tauSum > 0.0 && std::isfinite(tauSum)) {
            auto addCarrierFloor = [&](int row, Real carrierDensity, Real ni, Real volume, Real fallbackSign) {
                if (constrainedRows[static_cast<std::size_t>(row)] || ni <= 0.0 || volume <= 0.0)
                    return;
                if (!(carrierDensity < carrierDiagonalFloor_.minorityDensityRatio * ni))
                    return;
                const Real rawFloor = carrierDiagonalFloor_.scale * volume * ni *
                    sourceIntegralFactor / (tauSum * Vt_);
                const Real floor = scaleDerivative(row, row, rawFloor);
                if (!(floor > 0.0) || !std::isfinite(floor))
                    return;
                const Real diagonal = assembledDiagonal[static_cast<std::size_t>(row)];
                const Real absDiagonal = std::abs(diagonal);
                if (!(absDiagonal < floor))
                    return;
                const Real sign = diagonal < 0.0 ? -1.0 : (diagonal > 0.0 ? 1.0 : fallbackSign);
                const Real addition = sign * (floor - absDiagonal);
                triplets.emplace_back(row, row, addition);
                assembledDiagonal[static_cast<std::size_t>(row)] += addition;
            };

            for (Index i = 0; i < Nidx; ++i) {
                const int ii = static_cast<int>(i);
                const Real ni = ni_[i];
                // In a depleted SRH generation row, dR/dqF is O(ni / ((taun + taup) Vt));
                // multiplying by nodal volume gives an absolute Jacobian scale.
                addCarrierFloor(phinOffset() + ii, n(ii), ni, vol_[i], -1.0);
                addCarrierFloor(phipOffset() + ii, p(ii), ni, vol_[i], 1.0);
            }
        }
    }

    for (const auto& [node, value] : bcs.psi) {
        (void)value;
        triplets.emplace_back(psiOffset() + static_cast<int>(node),
                              psiOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phin) {
        (void)value;
        triplets.emplace_back(phinOffset() + static_cast<int>(node),
                              phinOffset() + static_cast<int>(node), 1.0);
    }
    for (const auto& [node, value] : bcs.phip) {
        (void)value;
        triplets.emplace_back(phipOffset() + static_cast<int>(node),
                              phipOffset() + static_cast<int>(node), 1.0);
    }

    SparseMatrixd J(3 * N, 3 * N);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

SparseMatrixd CoupledDDAssembler::finiteDifferenceJacobian(
    const VectorXd& x,
    const CoupledDDBoundaryConditions& bcs,
    Real relativeStep) const
{
    if (relativeStep <= 0.0)
        throw std::invalid_argument(
            "CoupledDDAssembler::finiteDifferenceJacobian: relativeStep must be positive.");

    const int M = x.size();
    const VectorXd r0 = residual(x, bcs);
    std::vector<Eigen::Triplet<double>> triplets;
    // Heuristic: for a 2-D mesh with M = 3*N unknowns and sparse connectivity
    // (average ~6-7 neighbours per node), the Jacobian has O(N * 7 * 9) ~= 63*N
    // non-zeros.  Reserving M * 7 avoids the M^2 allocation that would OOM for
    // any realistically-sized mesh while still keeping reallocations rare.
    triplets.reserve(static_cast<std::size_t>(M) * 7);

    // Minimum absolute perturbation to prevent h == 0 when x(col) == 0.
    constexpr Real minAbsStep = 1.0e-15;

    for (int col = 0; col < M; ++col) {
        const Real h = std::max(relativeStep * std::max(1.0, std::abs(x(col))), minAbsStep);
        VectorXd xp = x;
        xp(col) += h;
        const VectorXd rp = residual(xp, bcs);
        const VectorXd dr = (rp - r0) / h;
        for (int row = 0; row < M; ++row) {
            if (dr(row) != 0.0)
                triplets.emplace_back(row, col, dr(row));
        }
    }

    SparseMatrixd J(M, M);
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

std::string CoupledDDAssembler::impactIonizationConfigurationFingerprint() const
{
    return bvProcessConfigurationFingerprint(
        mobilityConfig_, impactIonizationConfig_);
}

std::string CoupledDDAssembler::impactIonizationActiveBranchFingerprint(
    const VectorXd& x) const
{
    ScopedPerformanceTimer timer("dd.active_branch_fingerprint");
    incrementPerformanceCounter("dd.active_branch_fingerprint_calls");
    const CoupledDDState state = unpack(x);
    const Real potentialScale = scaling_.enabled ? scaling_.V0 : 1.0;

    DDSolution solution;
    solution.psi = state.psi * potentialScale;
    solution.phin = state.phin * potentialScale;
    solution.phip = state.phip * potentialScale;
    solution.n = electronDensity(x);
    solution.p = holeDensity(x);

    const BVProcessProbeResult probe = evaluateBVProcessProbe(
        mesh_,
        doping_,
        solution,
        mobilityConfig_,
        impactIonizationConfig_,
        bandgapNarrowingConfig_,
        matdb_,
        Vt_ * constants::q / constants::kb,
        scaling_.enabled ? scaling_.fieldFromCoordinateDeltaFactor : 1.0);

    std::string canonical = probe.configurationFingerprint;
    for (const BVProcessProbeRecord& record : probe.records) {
        canonical.push_back('|');
        canonical += record.activeBranchFingerprint;
    }
    return fnv1a64(canonical);
}

} // namespace vela
