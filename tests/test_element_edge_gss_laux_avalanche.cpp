#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"

#include <array>
#include <vector>

using namespace vela;

namespace {

DeviceMesh makeSingleRightTriangle()
{
    DeviceMesh mesh;
    Node n0;
    n0.id = 0;
    n0.x = 0.0;
    n0.y = 0.0;
    mesh.addNode(n0);
    Node n1;
    n1.id = 1;
    n1.x = 1.0e-6;
    n1.y = 0.0;
    mesh.addNode(n1);
    Node n2;
    n2.id = 2;
    n2.x = 0.0;
    n2.y = 0.5e-6;
    mesh.addNode(n2);

    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);

    Region region;
    region.id = 0;
    region.name = "si";
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);

    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("Element-edge GSS Laux mode has an explicit canonical contract",
          "[impact][element_edge_gss_laux]")
{
    ImpactIonizationModelConfig config;
    config.model = "van_overstraeten";
    config.drivingForce = "quasi_fermi_gradient";
    config.generation = "current_density";
    config.currentApproximation = "element_edge_sg_gss_laux";
    config.quasiFermiGradientDiscretization = "cell_gradient";
    config.sourceMappingMode = "element_vertex_box_measure";

    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(config, "test"));
    REQUIRE(detail::usesElementEdgeGssLauxAvalancheSource(config));
    config.drivingForce = "electric_field";
    config.quasiFermiGradientDiscretization = "edge_difference";
    REQUIRE_NOTHROW(
        detail::validateImpactIonizationDrivingForce(config, "test"));

    config.sourceMappingMode = "node_F_node_alpha_node_G";
    REQUIRE_THROWS_WITH(
        detail::validateImpactIonizationDrivingForce(config, "test"),
        Catch::Matchers::ContainsSubstring(
            "element_edge_sg_gss_laux requires the canonical element-box configuration"));
}

TEST_CASE("GSS Laux element reconstruction keeps a zero-box diagonal current",
          "[impact][element_edge_gss_laux][geometry]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    const Cell& cell = mesh.getCell(0);
    const Point2 expected{3.0, -4.0};
    std::array<Real, 3> signedEdgeCurrent{};

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 =
            cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        const Point2 delta =
            detail::meshPoint(mesh, node1) - detail::meshPoint(mesh, node0);
        const Point2 tangent = delta / delta.norm();
        signedEdgeCurrent[static_cast<std::size_t>(localEdge)] =
            tangent.dot(expected);
    }

    const auto partialVolumes =
        detail::tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    REQUIRE(partialVolumes[1] == Catch::Approx(0.0).margin(1.0e-30));
    REQUIRE(signedEdgeCurrent[1] != Catch::Approx(0.0));

    const Point2 reconstructed = detail::gssLauxTri3CurrentVector(
        mesh, cell, signedEdgeCurrent);
    REQUIRE(reconstructed.x() ==
            Catch::Approx(expected.x()).epsilon(1.0e-13));
    REQUIRE(reconstructed.y() ==
            Catch::Approx(expected.y()).epsilon(1.0e-13));

    const auto vertexMeasures =
        detail::tri3ElementVertexBoxMeasures(mesh, cell);
    const Real measureSum =
        vertexMeasures[0] + vertexMeasures[1] + vertexMeasures[2];
    REQUIRE(measureSum == Catch::Approx(0.25e-12).epsilon(1.0e-13));
}

TEST_CASE("Element-edge GSS Laux records use exact SG currents and box measures",
          "[impact][element_edge_gss_laux][source]")
{
    DeviceMesh mesh = makeSingleRightTriangle();
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    MaterialDatabase matdb;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {RegionDopingSpec{"si", 1.0e21, 0.0}});
    const auto cellMaterials =
        detail::buildCellMaterials(mesh, matdb, constants::T0);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "element_edge_sg_gss_laux";
    impactConfig.quasiFermiGradientDiscretization = "cell_gradient";
    impactConfig.sourceMappingMode = "element_vertex_box_measure";
    const auto impact = makeImpactIonizationModel(impactConfig);

    MobilityModelConfig mobilityConfig = mobilityModelConfig("constant");
    const auto mobility = makeMobilityModel(mobilityConfig);
    VectorXd psi(mesh.numNodes());
    VectorXd phin(mesh.numNodes());
    VectorXd phip(mesh.numNodes());
    VectorXd n(mesh.numNodes());
    VectorXd p(mesh.numNodes());
    psi << 0.0, -0.20, 0.08;
    phin << 0.0, -0.40, 0.10;
    phip << 0.0, 0.30, -0.12;
    n << 1.0e20, 3.0e20, 1.5e20;
    p << 2.0e20, 6.0e20, 2.5e20;
    const std::vector<Real> ni(
        static_cast<std::size_t>(mesh.numNodes()), 1.0e16);
    const Real Vt = constants::kb * constants::T0 / constants::q;

    const auto record =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, *mobility, cellEdges.at(0), mesh, doping,
            cellMaterials, 0, psi, phin, phip, n, p, ni, Vt);

    REQUIRE(record.cellId == 0);
    REQUIRE(record.vertexMeasures[0] + record.vertexMeasures[1] +
                record.vertexMeasures[2] ==
            Catch::Approx(0.25e-12).epsilon(1.0e-13));
    REQUIRE(record.edgePartialVolumes[1] ==
            Catch::Approx(0.0).margin(1.0e-30));
    REQUIRE(record.electronSignedEdgeFlux[1] != Catch::Approx(0.0));
    REQUIRE(record.holeSignedEdgeFlux[1] != Catch::Approx(0.0));

    const Real totalMeasure = 0.25e-12;
    REQUIRE(record.electronCurrentVector.norm() > 0.0);
    REQUIRE(record.holeCurrentVector.norm() > 0.0);
    REQUIRE(record.electronSourceIntegrals[0] +
                record.electronSourceIntegrals[1] +
                record.electronSourceIntegrals[2] ==
            Catch::Approx(record.electronAlpha *
                          record.electronCurrentVector.norm() *
                          totalMeasure)
                .epsilon(1.0e-13));
    REQUIRE(record.holeSourceIntegrals[0] +
                record.holeSourceIntegrals[1] +
                record.holeSourceIntegrals[2] ==
            Catch::Approx(record.holeAlpha *
                          record.holeCurrentVector.norm() *
                          totalMeasure)
                .epsilon(1.0e-13));

    const auto components =
        detail::currentDensityAvalancheSourceComponentIntegrals(
            impactConfig, *impact, mobilityConfig, *mobility, edgeCells,
            mesh, doping, cellMaterials, psi, phin, phip, n, p, ni, Vt);
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        const Index node = mesh.getCell(0).node_ids[localNode];
        REQUIRE(components.electron[node] ==
                Catch::Approx(record.electronSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.hole[node] ==
                Catch::Approx(record.holeSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.combined[node] ==
                Catch::Approx(record.combinedSourceIntegrals[localNode])
                    .epsilon(1.0e-13));
        REQUIRE(components.electron[node] + components.hole[node] ==
                Catch::Approx(components.combined[node])
                    .epsilon(1.0e-13));
    }
    impactConfig.drivingForce = "electric_field";
    const auto electricRecord =
        detail::elementEdgeGssLauxAvalancheSourceRecordForCell(
            impactConfig, *impact, *mobility, cellEdges.at(0), mesh, doping,
            cellMaterials, 0, psi, phin, phip, n, p, ni, Vt);
    const Real expectedElectricField = std::hypot(2.0e5, 1.6e5);
    REQUIRE(electricRecord.electronImpactField ==
            Catch::Approx(expectedElectricField).epsilon(1.0e-13));
    REQUIRE(electricRecord.holeImpactField ==
            Catch::Approx(expectedElectricField).epsilon(1.0e-13));
    REQUIRE(electricRecord.electronImpactField !=
            Catch::Approx(record.electronImpactField));
}
