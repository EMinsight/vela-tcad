#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/BVProcessProbe.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"

#include <array>
#include <cmath>
#include <set>
#include <string>

using namespace vela;

namespace {

DeviceMesh makeTriangle(
    std::array<Point2, 3> points,
    bool reverseOrientation = false,
    bool addContact = false)
{
    DeviceMesh mesh;
    for (Index nodeId = 0; nodeId < 3; ++nodeId) {
        Node node;
        node.id = nodeId;
        node.x = points[nodeId].x();
        node.y = points[nodeId].y();
        mesh.addNode(node);
    }
    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = reverseOrientation
        ? std::vector<Index>{0, 2, 1}
        : std::vector<Index>{0, 1, 2};
    mesh.addCell(cell);
    Region region;
    region.id = 0;
    region.name = "silicon";
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);
    if (addContact) {
        Contact contact;
        contact.id = 0;
        contact.name = "anode";
        contact.region_id = 0;
        contact.node_ids = {0};
        mesh.addContact(contact);
    }
    mesh.buildEdges();
    return mesh;
}

DopingModel makeDoping(Index count)
{
    DopingModel doping(count);
    for (Index node = 0; node < count; ++node)
        doping.setNodeDoping(node, 1.0e21 + node * 1.0e20, 0.0);
    return doping;
}

DDSolution makeState(Index count)
{
    DDSolution state;
    state.psi.resize(static_cast<Eigen::Index>(count));
    state.phin.resize(static_cast<Eigen::Index>(count));
    state.phip.resize(static_cast<Eigen::Index>(count));
    state.n.resize(static_cast<Eigen::Index>(count));
    state.p.resize(static_cast<Eigen::Index>(count));
    for (Index node = 0; node < count; ++node) {
        const Eigen::Index i = static_cast<Eigen::Index>(node);
        state.psi(i) = 0.08 * static_cast<Real>(node);
        state.phin(i) = 0.03 * static_cast<Real>(node);
        state.phip(i) = -0.02 * static_cast<Real>(node);
        state.n(i) = 1.0e20 * static_cast<Real>(node + 1);
        state.p(i) = 1.0e20 * static_cast<Real>(count - node);
    }
    state.converged = true;
    return state;
}

ImpactIonizationModelConfig triangleImpact(std::string coupling)
{
    ImpactIonizationModelConfig config;
    config.model = "selberherr";
    config.couplingMode = std::move(coupling);
    config.drivingForce = "quasi_fermi_gradient";
    config.generation = "current_density";
    config.currentApproximation = "cell_reconstructed";
    config.quasiFermiGradientDiscretization = "cell_gradient";
    config.sourceMappingMode = "triangle_gss_gradqf_truncated";
    config.electronA = 1.0;
    config.electronB = 1.0;
    config.holeA = 1.0;
    config.holeB = 1.0;
    return config;
}

Real relativeDifference(Real first, Real second)
{
    return std::abs(first - second) /
        std::max({std::abs(first), std::abs(second), Real{1.0e-300}});
}

BVProcessProbeResult evaluate(
    const DeviceMesh& mesh,
    const DDSolution& state,
    const ImpactIonizationModelConfig& impact)
{
    return evaluateBVProcessProbe(
        mesh,
        makeDoping(mesh.numNodes()),
        state,
        mobilityModelConfig("constant"),
        impact,
        BandgapNarrowingConfig{},
        MaterialDatabase{},
        300.0);
}

} // namespace

TEST_CASE("BV process probe closes production source and residual scatter",
          "[bv_process_probe][closure]")
{
    const DeviceMesh mesh = makeTriangle(
        {Point2{0.0, 0.0}, Point2{1.0e-6, 0.0}, Point2{0.0, 1.0e-6}});
    DDSolution state = makeState(mesh.numNodes());
    const ImpactIonizationModelConfig impact =
        triangleImpact("self_consistent");
    CoupledDDAssembler assembler(
        mesh,
        MaterialDatabase{},
        makeDoping(mesh.numNodes()),
        constants::kb * 300.0 / constants::q,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impact);
    CoupledDDState coupledState;
    coupledState.psi = state.psi;
    coupledState.phin = state.phin;
    coupledState.phip = state.phip;
    const VectorXd x = assembler.pack(coupledState);
    state.n = assembler.electronDensity(x);
    state.p = assembler.holeDensity(x);
    const BVProcessProbeResult probe = evaluate(mesh, state, impact);

    REQUIRE(probe.records.size() == 6);
    REQUIRE(assembler.impactIonizationConfigurationFingerprint() ==
            probe.configurationFingerprint);
    REQUIRE(probe.totalSourceIntegral > 0.0);
    REQUIRE(relativeDifference(
                probe.totalQGContribution,
                constants::q * probe.totalSourceIntegral) <
            1.0e-15);
    REQUIRE(relativeDifference(
                probe.electronResidualContribution,
                -probe.totalSourceIntegral) <
            1.0e-15);
    REQUIRE(relativeDifference(
                probe.holeResidualContribution,
                -probe.totalSourceIntegral) <
            1.0e-15);

    const VectorXd sourceResidual = assembler.impactIonizationSourceResidual(
        x, CoupledDDBoundaryConditions{});
    Real electronResidual = 0.0;
    Real holeResidual = 0.0;
    const Eigen::Index nodeCount =
        static_cast<Eigen::Index>(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        electronResidual += sourceResidual(
            nodeCount + static_cast<Eigen::Index>(node));
        holeResidual += sourceResidual(
            2 * nodeCount + static_cast<Eigen::Index>(node));
    }
    REQUIRE(relativeDifference(
                electronResidual, probe.electronResidualContribution) <
            1.0e-12);
    REQUIRE(relativeDifference(
                holeResidual, probe.holeResidualContribution) <
            1.0e-12);

    for (const BVProcessProbeRecord& record : probe.records) {
        REQUIRE(record.configurationFingerprint ==
                probe.configurationFingerprint);
        REQUIRE(record.scatterCount == 2);
        REQUIRE(record.sourceWeights[0] + record.sourceWeights[1] ==
                Catch::Approx(1.0).epsilon(1.0e-15));
        REQUIRE_FALSE(record.activeBranchFingerprint.empty());
        REQUIRE(std::isfinite(record.mobilityLimiter));
        REQUIRE(std::isfinite(record.generationRate));
    }
}

TEST_CASE("BV process probe postprocess branch observes source without residual scatter",
          "[bv_process_probe][postprocess_only]")
{
    const DeviceMesh mesh = makeTriangle(
        {Point2{0.0, 0.0}, Point2{1.0e-6, 0.0}, Point2{0.2e-6, 0.7e-6}});
    const DDSolution state = makeState(mesh.numNodes());
    const BVProcessProbeResult coupled =
        evaluate(mesh, state, triangleImpact("self_consistent"));
    const BVProcessProbeResult observed =
        evaluate(mesh, state, triangleImpact("postprocess_only"));

    REQUIRE(observed.totalSourceIntegral ==
            Catch::Approx(coupled.totalSourceIntegral).epsilon(1.0e-15));
    REQUIRE(observed.totalQGContribution ==
            Catch::Approx(coupled.totalQGContribution).epsilon(1.0e-15));
    REQUIRE(observed.electronResidualContribution == 0.0);
    REQUIRE(observed.holeResidualContribution == 0.0);
    REQUIRE(observed.configurationFingerprint !=
            coupled.configurationFingerprint);
    for (const BVProcessProbeRecord& record : observed.records) {
        REQUIRE_FALSE(record.solverCoupled);
        REQUIRE(record.electronResidualContributions[0] == 0.0);
        REQUIRE(record.holeResidualContributions[0] == 0.0);
    }
}

TEST_CASE("BV process probe records geometry orientation contact and carrier branches",
          "[bv_process_probe][geometry]")
{
    const std::array<std::array<Point2, 3>, 3> shapes = {
        std::array<Point2, 3>{
            Point2{0.0, 0.0}, Point2{1.0e-6, 0.0}, Point2{0.5e-6, 0.8e-6}},
        std::array<Point2, 3>{
            Point2{0.0, 0.0}, Point2{1.0e-6, 0.0}, Point2{0.0, 1.0e-6}},
        std::array<Point2, 3>{
            Point2{0.0, 0.0}, Point2{1.0e-6, 0.0}, Point2{0.1e-6, 0.1e-6}},
    };
    for (const auto& shape : shapes) {
        for (bool reversed : {false, true}) {
            const DeviceMesh mesh = makeTriangle(shape, reversed, true);
            const BVProcessProbeResult probe =
                evaluate(mesh, makeState(3), triangleImpact("self_consistent"));
            REQUIRE(probe.records.size() == 6);
            std::set<std::string> carriers;
            bool sawContact = false;
            bool sawInteriorEdge = false;
            for (const BVProcessProbeRecord& record : probe.records) {
                carriers.insert(record.carrier);
                sawContact = sawContact || record.contactAdjacent;
                sawInteriorEdge = sawInteriorEdge || !record.contactAdjacent;
                REQUIRE(std::isfinite(record.sourceIntegral));
                REQUIRE(std::isfinite(record.currentVector.x()));
                REQUIRE(std::isfinite(record.currentVector.y()));
            }
            REQUIRE(carriers ==
                    std::set<std::string>{"electron", "hole"});
            REQUIRE(sawContact);
            REQUIRE(sawInteriorEdge);
        }
    }
}

TEST_CASE("BV process configuration fingerprint is deterministic and branch sensitive",
          "[bv_process_probe][fingerprint]")
{
    const MobilityModelConfig mobility = mobilityModelConfig("constant");
    const ImpactIonizationModelConfig self =
        triangleImpact("self_consistent");
    ImpactIonizationModelConfig post = self;
    post.couplingMode = "postprocess_only";
    REQUIRE(bvProcessConfigurationFingerprint(mobility, self) ==
            bvProcessConfigurationFingerprint(mobility, self));
    REQUIRE(bvProcessConfigurationFingerprint(mobility, self) !=
            bvProcessConfigurationFingerprint(mobility, post));
}
