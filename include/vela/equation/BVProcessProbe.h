#pragma once

#include "vela/core/Types.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/solver/GummelSolver.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace vela {

/// A normalized, solver-provenance record for one carrier contribution on one
/// production avalanche support. Values stay in the active solver unit system;
/// output adapters are responsible for display-unit conversion.
struct BVProcessProbeRecord {
    std::string supportKind;
    std::string carrier;
    Index cellId = 0;
    int localEdge = -1;
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;

    Real psi0 = 0.0;
    Real psi1 = 0.0;
    Real quasiFermi0 = 0.0;
    Real quasiFermi1 = 0.0;
    Real density0 = 0.0;
    Real density1 = 0.0;
    Real midpointDensity = 0.0;

    Point2 electricFieldVector = Point2::Zero();
    Point2 quasiFermiGradientVector = Point2::Zero();
    Point2 currentVector = Point2::Zero();
    std::string currentVectorProvenance;

    Real lowFieldMobility = 0.0;
    Real highFieldDrive = 0.0;
    Real finalMobility = 0.0;
    Real mobilityLimiter = 0.0;
    Real directedSgFlux = 0.0;
    Real selectedFluxMagnitude = 0.0;
    Real impactField = 0.0;
    Real alpha = 0.0;
    Real sourceMeasure = 0.0;
    Real generationRate = 0.0;
    Real sourceIntegral = 0.0;
    Real qGContribution = 0.0;

    std::array<Index, 6> scatterNodes{};
    std::array<Real, 6> sourceWeights{};
    std::array<Real, 6> electronResidualContributions{};
    std::array<Real, 6> holeResidualContributions{};
    std::size_t scatterCount = 0;

    bool solverCoupled = false;
    bool contactAdjacent = false;
    bool zeroMeasure = false;
    bool zeroMobility = false;
    bool zeroAlpha = false;
    bool reconstructedCurrent = false;
    bool directionalPartition = false;
    std::string productionBranchDetails;
    std::string activeBranches;
    std::string configurationFingerprint;
    std::string activeBranchFingerprint;
};

struct BVProcessProbeResult {
    std::vector<BVProcessProbeRecord> records;
    std::string configurationFingerprint;
    Real totalSourceIntegral = 0.0;
    Real totalQGContribution = 0.0;
    Real electronResidualContribution = 0.0;
    Real holeResidualContribution = 0.0;
};

std::string bvProcessConfigurationFingerprint(
    const MobilityModelConfig& mobility,
    const ImpactIonizationModelConfig& impact);

BVProcessProbeResult evaluateBVProcessProbe(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const DDSolution& state,
    const MobilityModelConfig& mobility,
    const ImpactIonizationModelConfig& impact,
    const BandgapNarrowingConfig& bandgapNarrowing,
    const MaterialDatabase& materials,
    Real temperature_K,
    Real fieldFromCoordinateDeltaFactor = 1.0);

} // namespace vela
