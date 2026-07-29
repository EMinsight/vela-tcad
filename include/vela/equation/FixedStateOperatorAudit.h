#pragma once

#include "vela/equation/BVProcessProbe.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/solver/GummelSolver.h"

#include <array>
#include <vector>

namespace vela {

using TriangleGssAvalancheSourceRecord = detail::TriangleGssAvalancheSourceRecord;
using ElementEdgeGssLauxAvalancheSourceRecord =
    detail::ElementEdgeGssLauxAvalancheSourceRecord;

struct FixedStateNodeRecord {
    Index nodeId;
    Real psi, phin, phip, n, p;
};

struct FixedStateEdgeRecord {
    Index edgeId, node0, node1;
    Real length, electronRawSignedFlux, holeRawSignedFlux;
    Real electronMidpointDensity, holeMidpointDensity;
    Real electronImpactField, holeImpactField;
    Real electronAlpha, holeAlpha, edgeArea;
};

struct FixedStateTriangleRecord {
    Index cellId;
    std::array<Index, 3> nodes;
    Real signedDoubleArea;
    Point2 gradPsi, gradPhin, gradPhip;
    std::vector<TriangleGssAvalancheSourceRecord> localEdges;
};

struct FixedStateOperatorAuditOptions {
    bool allowGeneralTri3 = false;
};

struct FixedStateOperatorAuditResult {
    std::vector<FixedStateNodeRecord> nodes;
    std::vector<FixedStateEdgeRecord> edges;
    std::vector<FixedStateTriangleRecord> triangles;
    std::vector<ElementEdgeGssLauxAvalancheSourceRecord>
        elementEdgeGssLauxTriangles;
    BVProcessProbeResult processProbe;
};

FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const VectorXd& doping,
    const DDSolution& state,
    const GummelConfig& config);

FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const VectorXd& doping,
    const DDSolution& state,
    const GummelConfig& config,
    const MaterialDatabase& materials);
FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const DDSolution& state,
    const GummelConfig& config,
    const MaterialDatabase& materials,
    FixedStateOperatorAuditOptions options = {});

} // namespace vela
