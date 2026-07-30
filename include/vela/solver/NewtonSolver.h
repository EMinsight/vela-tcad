#pragma once

#include "vela/core/PhysicalConstants.h"
#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/numerics/LineSearch.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/simulation/PseudoArclength.h"
#include "vela/solver/GummelSolver.h"
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace vela {

struct NewtonCarrierRowConvergenceConfig {
    std::string mode = "off"; ///< "off", "report", or "enforce".
    Real epsRow = 1.0e-3;
    Real scaleFloor = 1.0e-300;
    Real minSourceScaleFraction = 1.0e-3;
    Real minSourceScale = 0.0;
    int minEnforceMaxIter = 200;
    std::string diagnosticCsvFile;
    std::string traceCsvFile;
    std::vector<Index> traceNodes;
    int traceFirstIterations = 10;
    int traceEveryIterations = 10;
};

struct NewtonContinuityRowScalingConfig {
    bool enabled = false;
    Real fluxFraction = 1.0e-3;
    Real scaleFloor = 1.0e-30;
    Real minSourceScale = 0.0;
    Real minWeight = 1.0e-12;
    Real maxWeight = 1.0e12;
};

struct NewtonGlobalContinuityClosureConfig {
    std::string mode = "off"; ///< "off", "report", or "enforce".
    Real tolerance = 1.0e-2;
    Real sourceFloor = 1.0e-12;
};

struct NewtonGlobalContinuityCarrierClosure {
    bool qualified = false;
    Real contactFlux = 0.0;
    Real integratedSource = 0.0;
    Real mismatch = 0.0;
    Real ratio = 0.0;
};

struct NewtonGlobalContinuityClosureEvaluation {
    bool enabled = false;
    bool enforced = false;
    bool satisfied = true;
    NewtonGlobalContinuityCarrierClosure electron;
    NewtonGlobalContinuityCarrierClosure hole;
};

struct NewtonCarrierRowConvergenceViolation {
    Index nodeId = 0;
    std::string carrier;
    Real residual = 0.0;
    Real scale = 0.0;
    Real ratio = 0.0;
    Real flux = 0.0;
    Real recombination = 0.0;
    Real impact = 0.0;
};

struct NewtonCarrierRowConvergenceEvaluation {
    bool enabled = false;
    bool enforced = false;
    bool satisfied = true;
    Real epsRow = 0.0;
    Real maxRatio = 0.0;
    Index maxRatioNode = -1;
    std::string maxRatioCarrier;
    std::vector<NewtonCarrierRowConvergenceViolation> violations;
};

struct NewtonCarrierRowRecoveryConfig {
    std::string mode = "off"; ///< "off" or "gummel_density".
    int maxAttempts = 1; ///< Maximum density passes per recovery cycle.
    int maxCycles = 1; ///< Maximum Newton/recovery cycles at one bias point.
    Real densityChangeReltol = 1.0e-8; ///< Stop density passes when max relative density change falls below this value.
};

struct NewtonCarrierRowRecoveryResult {
    DDSolution solution;
    bool attempted = false;
    std::string mode;
    int electronRowsUpdated = 0;
    int holeRowsUpdated = 0;
    int densityPasses = 0;
    int cyclesAttempted = 0;
    bool densityConverged = false;
    Real maxDensityRelativeChange = 0.0;
    Real maxPsiDelta_V = 0.0;
    Real maxCarrierDensityRatio = 0.0;
};
struct NewtonConfig {
    int maxIter = 20;
    Real reltol = 1.0e-8;
    Real abstol = 1.0e-18;
    Real temperature_K = constants::T0; ///< Lattice temperature [K]
    Real dampingFactor = 1.0;
    bool lineSearch = true;
    bool verbose = true;
    bool warmStart = false; ///< Preserve supplied quasi-Fermi potentials instead of resetting interiors.
    bool diagnostics = false; ///< Store detailed line-search diagnostics in NewtonResult::history.
    Real maxUpdate = 0.0; ///< Optional infinity-norm cap on one Newton update in solver unknown units.
    Real quasiFermiUpdateLimit_V = 0.0; ///< Optional physical-voltage cap on phin/phip Newton updates.
    Real quasiFermiUpdateLimitMinority_V = 0.0; ///< Optional tighter physical-voltage cap on the minority-carrier quasi-Fermi update per node; 0 disables (uses the global cap for both carriers).
    Real stallResidualFloor = 1.0e-9; ///< Residual floor for accepting line-search stalls as solved.
    Real poissonLineSearchStallResidualFloor = 1.0e-6; ///< Poisson-block floor for near-flat line-search stalls.
    Real poissonLineSearchStallRelativeIncrease = 1.0e-5; ///< Allowed best rejected residual increase at the Poisson floor.
    Real poissonLineSearchStallCarrierResidualFloor = 1.0e-6; ///< Carrier-block ceiling for Poisson-floor stall acceptance.
    Real poissonLineSearchStallContactMajorityQfDropLimit_V = 5.0e-11; ///< Maximum contact-edge majority-carrier quasi-Fermi drop allowed for Poisson-floor stall acceptance; 0 disables.
    Real carrierRegularizationScale = 0.0; ///< Optional carrier-row diagonal regularization scale.
    CarrierDiagonalFloorRegularizationConfig carrierDiagonalFloor{}; ///< Optional absolute floor for depleted minority carrier-row diagonals.
    NewtonCarrierRowConvergenceConfig carrierRowConvergence{}; ///< Optional per-carrier-row local residual convergence check.
    NewtonCarrierRowRecoveryConfig carrierRowRecovery{}; ///< Optional recovery pass for locally unbalanced carrier rows.
    NewtonContinuityRowScalingConfig continuityRowScaling{}; ///< Optional source-aware left row equilibration.
    NewtonGlobalContinuityClosureConfig globalContinuityClosure{}; ///< Optional contact-flux versus integrated-source convergence check.
    Real finiteDifferenceStep = 1.0e-6;
    std::string jacobian = "analytic"; ///< "analytic" or "finite_difference"
    std::string quasiFermiReference = "none"; ///< "none" or "contact_majority"
    std::string residualNorm = "block"; ///< "block" or "l2" convergence/line-search norm
    Real residualWeightPsi = 1.0;
    Real residualWeightPhin = 1.0;
    Real residualWeightPhip = 1.0;
    Real residualScalePsi = 0.0;  ///< <= 0 selects max(initial psi-block residual norm, 1)
    Real residualScalePhin = 0.0; ///< <= 0 selects max(initial electron-continuity residual norm, 1)
    Real residualScalePhip = 0.0; ///< <= 0 selects max(initial hole-continuity residual norm, 1)
    std::string contactBoundaryReconstruction = "dominant_signed_contact_mean";
    bool contactBoundaryMinorityElectronRelaxation = true;
    Real contactBoundaryMinorityElectronRelaxationBiasThreshold_V = 0.1;
    bool contactBoundaryMinorityElectronRelaxationTwoTerminalOnly = true;
    std::string contactBoundaryMinorityElectronRelaxationContactSide = "p_contact_only";
    Real contactBoundaryMinorityElectronRelaxationStrength = 1.0;
    UnitScalingConfig inputScaling{}; ///< Input-unit mode from top-level config.
    UnitScalingReferenceConfig unitScalingRefs{}; ///< Optional reference overrides.
    Real taun = 1.0e-5;
    Real taup = 3.0e-6;
    Real augerCn = 2.90e-43; ///< Electron Auger coefficient [m^6/s]
    Real augerCp = 1.028e-43; ///< Hole Auger coefficient [m^6/s]
    MobilityModelConfig mobility{}; ///< Mobility model configuration
    std::vector<std::string> recombination = {"srh"}; ///< e.g. {"srh", "auger"}
    ImpactIonizationModelConfig impactIonization; ///< Avalanche generation model.
    BandgapNarrowingConfig bandgapNarrowing; ///< Effective ni model for high doping.
};

struct NewtonBlockResidualInfo {
    Real psi = 0.0;
    Real phin = 0.0;
    Real phip = 0.0;
    Real combined = 0.0;
};

struct NewtonResidualPeak {
    Index nodeId = 0;
    Real signedResidual = 0.0;
    Real absoluteResidual = 0.0;
};

struct NewtonIterationInfo {
    int iter = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
    Real relativeResidualNorm = 0.0;
    Real rawStepNorm = 0.0;
    int lineSearchAttempts = 0;
    bool lineSearchAccepted = false;
    NewtonBlockResidualInfo blockResiduals;
    NewtonBlockResidualInfo rowScaledBlockResiduals;
    NewtonResidualPeak topPoissonResidual;
    NewtonResidualPeak topElectronResidual;
    NewtonResidualPeak topHoleResidual;
    std::string sourceJacobianActiveBranchFingerprint;
    std::string event;
    NewtonCarrierRowConvergenceEvaluation carrierRowConvergence;
    std::vector<LineSearchIterationInfo> lineSearchHistory;
};

struct NewtonCarrierDiagnostics {
    bool positiveFinite = true;
    Real minElectronDensity = 0.0;
    Real minHoleDensity = 0.0;
    int nonfiniteElectronCount = 0;
    int nonfiniteHoleCount = 0;
    int nonpositiveElectronCount = 0;
    int nonpositiveHoleCount = 0;
};

struct NewtonTopResidualNode {
    Index nodeId = 0;
    Real x = 0.0;
    Real y = 0.0;
    Real poissonResidual = 0.0;
    Real absPoissonResidual = 0.0;
    Real donors = 0.0;
    Real acceptors = 0.0;
    Real netDoping = 0.0;
    Real effectiveIntrinsicDensity = 0.0;
};

struct NewtonFailureDiagnostics {
    std::string failureReason;
    int failedIteration = 0;
    Real residualNorm = 0.0;
    Real stepNorm = 0.0;
    Real dampingFactor = 0.0;
    int lineSearchAttempts = 0;
    std::string lineSearchFailureReason;
    NewtonBlockResidualInfo blockResiduals;
    NewtonCarrierDiagnostics carrierDiagnostics;
    Real maxContactMajorityQfDrop = 0.0;
    Real bestRejectedContactMajorityQfDrop = 0.0;
    std::vector<LineSearchIterationInfo> lineSearchHistory;
    std::vector<NewtonTopResidualNode> topPoissonResidualNodes;
};

struct NewtonResult {
    DDSolution solution;
    bool converged = false;
    int iters = 0;
    Real initialResidualNorm = 0.0;
    Real finalResidualNorm = 0.0;
    std::string convergenceReason;
    NewtonBlockResidualInfo finalBlockNorms;
    NewtonCarrierRowConvergenceEvaluation finalCarrierRowConvergence;
    NewtonGlobalContinuityClosureEvaluation finalGlobalContinuityClosure;
    NewtonCarrierRowRecoveryResult carrierRowRecovery;
    std::vector<NewtonIterationInfo> history;
    std::vector<NewtonIterationInfo> trace;
    NewtonFailureDiagnostics failureDiagnostics;
};

struct NewtonResidualEvaluation {
    VectorXd raw;
    NewtonBlockResidualInfo blockNorms;
    std::vector<Real> intrinsicDensity;
    bool scaledState = false;
    Real potentialScale = 1.0;
};

struct NewtonStepEvaluation {
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonDirectionalDerivativeEvaluation {
    NewtonResidualEvaluation residual;
    VectorXd perturbationPsi;
    VectorXd perturbationPhin;
    VectorXd perturbationPhip;
    VectorXd analyticJv;
    VectorXd finiteDifferenceJv;
    VectorXd forwardResidual;
    VectorXd backwardResidual;
    Real perturbationNorm = 0.0;
    Real analyticNorm = 0.0;
    Real finiteDifferenceNorm = 0.0;
    Real absoluteError = 0.0;
    Real relativeError = 0.0;
};

struct NewtonBlockStepEvaluation {
    std::string mode;
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonPoissonBlockInitialization {
    DDSolution coldInitial;
    DDSolution poissonBlockInitial;
    NewtonBlockResidualInfo coldBlockResiduals;
    NewtonBlockResidualInfo poissonBlockResiduals;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};

struct NewtonRegularizedCarrierStepEvaluation {
    Real regularizationScale = 0.0;
    NewtonResidualEvaluation residual;
    NewtonResidualEvaluation trialResidual;
    DDSolution trialSolution;
    VectorXd deltaPsi;
    VectorXd deltaPhin;
    VectorXd deltaPhip;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
    Real regularizationDiagonalNorm = 0.0;
};

struct NewtonCarrierRowDiagnostic {
    Index nodeId = 0;
    Real electronResidual = 0.0;
    Real holeResidual = 0.0;
    Real electronDiagonal = 0.0;
    Real holeDiagonal = 0.0;
    Real electronRowAbsSum = 0.0;
    Real holeRowAbsSum = 0.0;
    Real electronOffdiagAbsSum = 0.0;
    Real holeOffdiagAbsSum = 0.0;
    Real electronRowL2Norm = 0.0;
    Real holeRowL2Norm = 0.0;
    Real rawDeltaPhin_V = 0.0;
    Real rawDeltaPhip_V = 0.0;
    Real cappedDeltaPhin_V = 0.0;
    Real cappedDeltaPhip_V = 0.0;
};

struct NewtonCarrierRowDiagnosticsEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<NewtonCarrierRowDiagnostic> rows;
    Real potentialScale = 1.0;
    Real rawCarrierStepNorm = 0.0;
    Real cappedCarrierStepNorm = 0.0;
};

struct NewtonCarrierTermDiagnosticsEvaluation {
    NewtonResidualEvaluation residual;
    std::vector<CoupledDDCarrierTermDiagnostic> rows;
};

struct NewtonJacobianBlockAuditRow {
    std::string block;
    std::string configurationFingerprint;
    std::string activeBranchFingerprint;
    Real analyticNorm = 0.0;
    Real fdNorm = 0.0;
    Real diffNorm = 0.0;
    Real relDiff = 0.0;
    Real analyticPsiColumnNorm = 0.0;
    Real fdPsiColumnNorm = 0.0;
    Real diffPsiColumnNorm = 0.0;
    Real relPsiColumnDiff = 0.0;
    Real analyticPhinColumnNorm = 0.0;
    Real fdPhinColumnNorm = 0.0;
    Real diffPhinColumnNorm = 0.0;
    Real relPhinColumnDiff = 0.0;
    Real analyticPhipColumnNorm = 0.0;
    Real fdPhipColumnNorm = 0.0;
    Real diffPhipColumnNorm = 0.0;
    Real relPhipColumnDiff = 0.0;
};

class NewtonSolver {
public:
    NewtonSolver(const DeviceMesh& mesh,
                 const MaterialDatabase& matdb,
                 const DopingModel& doping,
                 const std::unordered_map<std::string, Real>& contactBiases,
                 NewtonConfig cfg = {},
                 std::vector<RegionFixedChargeSpec> fixedCharges = {},
                 std::vector<InterfaceSheetChargeSpec> sheetCharges = {});

    NewtonResult solve() const;
    NewtonPoissonBlockInitialization buildPoissonBlockInitialization() const;
    NewtonResult solve(const DDSolution& initial) const;
    NewtonResidualEvaluation evaluateResidual(const DDSolution& state) const;
    Real maxContactMajorityQuasiFermiDrop(const DDSolution& state) const;
    NewtonStepEvaluation evaluateStep(const DDSolution& state) const;
    NewtonDirectionalDerivativeEvaluation evaluateDirectionalDerivative(
        const DDSolution& state,
        const DDSolution& physicalPerturbation) const;
    NewtonBlockStepEvaluation evaluateBlockStep(
        const DDSolution& state,
        const std::string& mode) const;
    NewtonRegularizedCarrierStepEvaluation evaluateRegularizedCarrierStep(
        const DDSolution& state,
        Real regularizationScale) const;
    NewtonCarrierRowDiagnosticsEvaluation evaluateCarrierRowDiagnostics(
        const DDSolution& state) const;
    NewtonCarrierTermDiagnosticsEvaluation evaluateCarrierTermDiagnostics(
        const DDSolution& state) const;
    std::vector<NewtonJacobianBlockAuditRow> evaluateJacobianBlockAudit(
        const DDSolution& state,
        Real finiteDifferenceStep = 1.0e-7,
        std::vector<std::string> blocks = {},
        const std::string& finiteDifferenceMode = "double_symmetric") const;
    std::vector<CoupledDDEdgeFluxDiagnostic> evaluateSgEdgeFluxDiagnostics(
        const DDSolution& state) const;

    /// Build a pseudo-arclength continuation system over the coupled drift-diffusion
    /// residual, using the bias on `activeContact` (in volts) as the continuation
    /// parameter lambda. The returned callbacks operate on the scaled packed state
    /// vector [psi, phin, phip] produced by CoupledDDAssembler::pack, reusing the
    /// exact assembler, boundary-condition construction, and Jacobian assembly of
    /// the standard Newton solve. `biasFiniteDifferenceStep_V` sizes the central
    /// finite difference used to estimate dF/dV.
    ArclengthSystem makeArclengthSystem(const std::string& activeContact,
                                        Real biasFiniteDifferenceStep_V = 1.0e-4) const;

    /// Convert a physical-unit DDSolution into the scaled packed state vector and
    /// back, consistent with makeArclengthSystem. These help drive the continuation
    /// from a converged Newton solution and interpret the corrected state.
    VectorXd packArclengthState(const DDSolution& state) const;
    DDSolution unpackArclengthState(const VectorXd& x) const;

private:
    void configureQuasiFermiReferences(CoupledDDAssembler& assembler) const;
    CoupledDDBoundaryConditions buildBoundaryConditions(
        const CoupledDDAssembler& assembler) const;
    CoupledDDBoundaryConditions buildBoundaryConditions(
        const CoupledDDAssembler& assembler,
        const std::unordered_map<std::string, Real>& contactBiases) const;
    std::shared_ptr<CoupledDDAssembler> makeArclengthAssembler() const;
    DDScalingSpec buildScalingSpec() const;
    DDSolution buildInitialGuess(const CoupledDDAssembler& assembler,
                                 const CoupledDDBoundaryConditions& bcs) const;
    DDSolution makeSolution(const CoupledDDAssembler& assembler,
                            const VectorXd& x,
                            int iters) const;

    const DeviceMesh& mesh_;
    const MaterialDatabase& matdb_;
    const DopingModel& doping_;
    std::unordered_map<std::string, Real> contactBiases_;
    NewtonConfig cfg_;
    std::vector<RegionFixedChargeSpec> fixedCharges_;
    std::vector<InterfaceSheetChargeSpec> sheetCharges_;
};

NewtonConfig newtonConfigFromJson(
    const nlohmann::json& cfg,
    UnitScalingConfig scaling = {});

NewtonCarrierRowConvergenceEvaluation evaluateCarrierRowConvergence(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const NewtonCarrierRowConvergenceConfig& cfg);

NewtonGlobalContinuityClosureEvaluation evaluateGlobalContinuityClosure(
    const std::vector<CoupledDDCarrierTermDiagnostic>& rows,
    const std::vector<Index>& electronContactNodes,
    const std::vector<Index>& holeContactNodes,
    const NewtonGlobalContinuityClosureConfig& cfg);

NewtonCarrierRowRecoveryResult recoverCarrierRowsWithGummelDensity(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const std::unordered_map<std::string, Real>& contactBiases,
    const NewtonConfig& cfg,
    const DDSolution& state,
    const std::vector<NewtonCarrierRowConvergenceViolation>& violations,
    const NewtonCarrierRowRecoveryConfig& recovery);

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg = {});

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg = {});


NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

NewtonResult runNewton(const DeviceMesh& mesh,
                       const MaterialDatabase& matdb,
                       const DopingModel& doping,
                       const std::unordered_map<std::string, Real>& contactBiases,
                       const DDSolution& initial,
                       const NewtonConfig& cfg,
                       std::vector<RegionFixedChargeSpec> fixedCharges,
                       std::vector<InterfaceSheetChargeSpec> sheetCharges);

} // namespace vela
