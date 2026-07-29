#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/equation/FixedStateOperatorAudit.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/MobilityModel.h"

#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

using namespace vela;

namespace {

DeviceMesh makeMinimal6Mesh(bool reverseFirstTriangle = false,
                            bool unexpectedConnectivity = false)
{
    DeviceMesh mesh;
    const std::array<Point2, 6> points = {
        Point2{0.0, 0.5e-6}, Point2{1.0e-6, 0.5e-6}, Point2{2.0e-6, 0.5e-6},
        Point2{2.0e-6, 0.0}, Point2{0.0, 0.0}, Point2{1.0e-6, 0.0},
    };
    for (Index nodeId = 0; nodeId < points.size(); ++nodeId) {
        Node node;
        node.id = nodeId;
        node.x = points[nodeId].x();
        node.y = points[nodeId].y();
        mesh.addNode(node);
    }

    const std::array<std::array<Index, 3>, 4> triangles = {
        std::array<Index, 3>{0, 4, 1},
        std::array<Index, 3>{4, 5, 1},
        std::array<Index, 3>{1, 5, 3},
        std::array<Index, 3>{1, 3, 2},
    };
    for (Index cellId = 0; cellId < triangles.size(); ++cellId) {
        Cell cell;
        cell.id = cellId;
        cell.type = CellType::Tri3;
        cell.region_id = 0;
        cell.node_ids.assign(triangles[cellId].begin(), triangles[cellId].end());
        if (reverseFirstTriangle && cellId == 0)
            std::swap(cell.node_ids[1], cell.node_ids[2]);
        if (unexpectedConnectivity && cellId == 3)
            cell.node_ids = {5, 3, 2};
        mesh.addCell(cell);
    }

    Region region;
    region.id = 0;
    region.name = "silicon";
    region.material = "Si";
    region.cell_ids = {0, 1, 2, 3};
    mesh.addRegion(region);
    mesh.buildEdges();
    return mesh;
}

DDSolution makeState()
{
    DDSolution state;
    state.psi.resize(6);
    state.phin.resize(6);
    state.phip.resize(6);
    state.n.resize(6);
    state.p.resize(6);
    state.psi << -0.20, -0.10, 0.00, -0.15, -0.05, 0.05;
    state.phin << -0.18, -0.08, 0.02, -0.13, -0.03, 0.07;
    state.phip << -0.22, -0.12, -0.02, -0.17, -0.07, 0.03;
    state.n << 1.0e20, 2.0e20, 3.0e20, 1.5e20, 2.5e20, 3.5e20;
    state.p << 3.5e20, 2.5e20, 1.5e20, 3.0e20, 2.0e20, 1.0e20;
    state.iters = 17;
    state.converged = false;
    return state;
}

VectorXd makeDoping()
{
    VectorXd doping(6);
    doping << -1.0e22, -5.0e21, 0.0, 0.0, 5.0e21, 1.0e22;
    return doping;
}

GummelConfig makeConfig()
{
    GummelConfig config;
    config.mobility = mobilityModelConfig("constant");
    config.impactIonization.model = "van_overstraeten";
    config.impactIonization.drivingForce = "quasi_fermi_gradient";
    config.impactIonization.generation = "current_density";
    config.impactIonization.currentApproximation = "density_gradient";
    config.impactIonization.quasiFermiGradientDiscretization = "cell_gradient";
    return config;
}

std::vector<unsigned char> vectorBytes(const VectorXd& value)
{
    const auto* first = reinterpret_cast<const unsigned char*>(value.data());
    return {first, first + static_cast<std::size_t>(value.size()) * sizeof(Real)};
}

} // namespace

TEST_CASE("fixed-state audit preserves supplied carrier state", "[minimal6][fixed-state]")
{
    const DeviceMesh mesh = makeMinimal6Mesh();
    const VectorXd doping = makeDoping();
    DDSolution state = makeState();
    const auto psiBytes = vectorBytes(state.psi);
    const auto phinBytes = vectorBytes(state.phin);
    const auto phipBytes = vectorBytes(state.phip);
    const auto nBytes = vectorBytes(state.n);
    const auto pBytes = vectorBytes(state.p);

    const auto result = evaluateFixedStateOperators(mesh, doping, state, makeConfig());

    REQUIRE(result.nodes.size() == 6);
    REQUIRE(result.nodes[1].psi == Catch::Approx(state.psi[1]).margin(0.0));
    REQUIRE(result.nodes[1].n == Catch::Approx(state.n[1]).margin(0.0));
    REQUIRE(vectorBytes(state.psi) == psiBytes);
    REQUIRE(vectorBytes(state.phin) == phinBytes);
    REQUIRE(vectorBytes(state.phip) == phipBytes);
    REQUIRE(vectorBytes(state.n) == nBytes);
    REQUIRE(vectorBytes(state.p) == pBytes);
    REQUIRE(state.iters == 17);
    REQUIRE_FALSE(state.converged);
}

TEST_CASE("fixed-state audit enumerates nine edges and four triangles",
          "[minimal6][fixed-state]")
{
    const DeviceMesh mesh = makeMinimal6Mesh();
    const auto result = evaluateFixedStateOperators(
        mesh, makeDoping(), makeState(), makeConfig());

    REQUIRE(result.edges.size() == 9);
    REQUIRE(result.triangles.size() == 4);
    REQUIRE(result.processProbe.records.size() == 14);
    REQUIRE(result.processProbe.totalSourceIntegral >= 0.0);
    REQUIRE(result.processProbe.electronResidualContribution ==
            Catch::Approx(-result.processProbe.totalSourceIntegral)
                .epsilon(1.0e-14));
    std::vector<Index> zeroCoupleEdges;
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        if (mesh.getEdge(edgeId).couple <= 0.0)
            zeroCoupleEdges.push_back(edgeId);
    }
    REQUIRE(mesh.numEdges() - zeroCoupleEdges.size() == 7);
    REQUIRE(zeroCoupleEdges == std::vector<Index>{1, 6});
    for (Index edgeId : zeroCoupleEdges) {
        REQUIRE(result.edges[edgeId].edgeId == edgeId);
        REQUIRE(std::isfinite(result.edges[edgeId].electronRawSignedFlux));
        REQUIRE(std::isfinite(result.edges[edgeId].holeRawSignedFlux));
    }
    REQUIRE(result.triangles[0].signedDoubleArea > 0.0);
    REQUIRE(result.triangles[0].nodes == std::array<Index, 3>{0, 4, 1});
    for (Index edgeId = 0; edgeId < result.edges.size(); ++edgeId) {
        REQUIRE(result.edges[edgeId].edgeId == edgeId);
        REQUIRE(result.edges[edgeId].node0 < result.edges[edgeId].node1);
        REQUIRE(std::isfinite(result.edges[edgeId].electronRawSignedFlux));
        REQUIRE(std::isfinite(result.edges[edgeId].holeRawSignedFlux));
    }
}

TEST_CASE("fixed-state audit honors an explicit material database",
          "[minimal6][fixed-state][materials]")
{
    const DeviceMesh mesh = makeMinimal6Mesh();
    const VectorXd doping = VectorXd::Zero(6);
    const DDSolution state = makeState();
    GummelConfig config = makeConfig();
    config.impactIonization.model = "none";

    MaterialDatabase baseline;
    MaterialDatabase overridden;
    Material silicon = overridden.getMaterial("Si");
    silicon.ni *= 2.0;
    overridden.addMaterial(silicon);

    const auto baselineResult =
        evaluateFixedStateOperators(mesh, doping, state, config, baseline);
    const auto overriddenResult =
        evaluateFixedStateOperators(mesh, doping, state, config, overridden);

    REQUIRE(baselineResult.edges[0].electronRawSignedFlux != 0.0);
    REQUIRE(overriddenResult.edges[0].electronRawSignedFlux ==
            Catch::Approx(2.0 * baselineResult.edges[0].electronRawSignedFlux)
                .epsilon(2.0e-14));
    REQUIRE(overriddenResult.edges[0].holeRawSignedFlux ==
            Catch::Approx(2.0 * baselineResult.edges[0].holeRawSignedFlux)
                .epsilon(2.0e-14));
}

TEST_CASE("fixed-state audit zeroes right-triangle hypotenuse support",
          "[minimal6][fixed-state][geometry]")
{
    const DeviceMesh mesh = makeMinimal6Mesh();
    const auto result = evaluateFixedStateOperators(
        mesh, makeDoping(), makeState(), makeConfig());

    REQUIRE(result.triangles[3].nodes == std::array<Index, 3>{1, 3, 2});
    REQUIRE(result.triangles[3].localEdges[0].node0 == 1);
    REQUIRE(result.triangles[3].localEdges[0].node1 == 3);
    REQUIRE(result.triangles[3].localEdges[0].truncatedPartialVolume == 0.0);
}

TEST_CASE("fixed-state audit exposes opt-in element-edge GSS Laux records",
          "[minimal6][fixed-state][element-edge]")
{
    GummelConfig config = makeConfig();
    config.impactIonization.currentApproximation =
        "element_edge_sg_gss_laux";
    config.impactIonization.sourceMappingMode =
        "element_vertex_box_measure";
    const auto result = evaluateFixedStateOperators(
        makeMinimal6Mesh(), makeDoping(), makeState(), config);

    REQUIRE(result.elementEdgeGssLauxTriangles.size() == 4);
    REQUIRE(result.processProbe.records.size() == 48);
    std::set<std::string> processSupports;
    for (const auto& process : result.processProbe.records)
        processSupports.insert(process.supportKind);
    REQUIRE(processSupports == std::set<std::string>{
        "element_edge_gss_laux", "element_vertex_gss_laux"});
    for (const auto& record : result.elementEdgeGssLauxTriangles) {
        REQUIRE(record.cellId < 4);
        REQUIRE(record.electronCurrentVector.allFinite());
        REQUIRE(record.holeCurrentVector.allFinite());
        REQUIRE(record.vertexMeasures[0] + record.vertexMeasures[1] +
                    record.vertexMeasures[2] ==
                Catch::Approx(0.25e-12).epsilon(1.0e-13));
        REQUIRE(record.electronSourceIntegrals[0] >= 0.0);
        REQUIRE(record.holeSourceIntegrals[0] >= 0.0);
    }
}

TEST_CASE("fixed-state audit explicitly accepts general Tri3 connectivity",
          "[fixed-state][general-tri3]")
{
    const DeviceMesh mesh = makeMinimal6Mesh(false, true);
    const VectorXd netDoping = makeDoping();
    DopingModel doping(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const Real value = netDoping(static_cast<Eigen::Index>(node));
        doping.setNodeDoping(node, std::max(value, 0.0), std::max(-value, 0.0));
    }
    MaterialDatabase materials;
    const auto result = evaluateFixedStateOperators(
        mesh, doping, makeState(), makeConfig(), materials,
        FixedStateOperatorAuditOptions{true});
    REQUIRE(result.nodes.size() == mesh.numNodes());
    REQUIRE(result.edges.size() == mesh.numEdges());
    REQUIRE(result.triangles.size() == mesh.numCells());
}
TEST_CASE("fixed-state audit rejects invalid contracts", "[minimal6][fixed-state]")
{
    const DeviceMesh mesh = makeMinimal6Mesh();
    const VectorXd doping = makeDoping();
    const GummelConfig config = makeConfig();

    SECTION("reversed triangle orientation") {
        REQUIRE_THROWS_WITH(
            evaluateFixedStateOperators(
                makeMinimal6Mesh(true), doping, makeState(), config),
            Catch::Matchers::ContainsSubstring("counter-clockwise"));
    }

    SECTION("unexpected connectivity") {
        REQUIRE_THROWS_WITH(
            evaluateFixedStateOperators(
                makeMinimal6Mesh(false, true), doping, makeState(), config),
            Catch::Matchers::ContainsSubstring("unexpected minimal6 connectivity"));
    }

    SECTION("wrong state vector size") {
        DDSolution state = makeState();
        state.p.conservativeResize(5);
        REQUIRE_THROWS_WITH(
            evaluateFixedStateOperators(mesh, doping, state, config),
            Catch::Matchers::ContainsSubstring("state vector size"));
    }

    SECTION("non-finite state") {
        DDSolution state = makeState();
        state.phin[2] = std::numeric_limits<Real>::infinity();
        REQUIRE_THROWS_WITH(
            evaluateFixedStateOperators(mesh, doping, state, config),
            Catch::Matchers::ContainsSubstring("finite"));
    }

    SECTION("wrong doping size") {
        REQUIRE_THROWS_WITH(
            evaluateFixedStateOperators(mesh, VectorXd::Zero(5), makeState(), config),
            Catch::Matchers::ContainsSubstring("doping size"));
    }
}
