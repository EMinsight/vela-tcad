#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using namespace vela;

namespace {

nlohmann::json readJson(const std::filesystem::path& path)
{
    std::ifstream input(path);
    REQUIRE(input.is_open());
    nlohmann::json json;
    input >> json;
    return json;
}

DopingModel dopingFromDeck(const DeviceMesh& mesh, const nlohmann::json& cfg)
{
    std::vector<RegionDopingSpec> specs;
    for (const auto& entry : cfg.at("doping")) {
        specs.push_back({
            entry.at("region").get<std::string>(),
            entry.at("donors").get<Real>(),
            entry.at("acceptors").get<Real>(),
        });
    }
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

bool vectorIsFinite(const VectorXd& values)
{
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values(i)))
            return false;
    }
    return true;
}

std::vector<std::vector<Index>> buildNodeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeCells(mesh.numNodes());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        for (Index nodeId : mesh.getCell(cellId).node_ids)
            nodeCells.at(nodeId).push_back(cellId);
    }
    return nodeCells;
}

std::vector<Index> contactNodesOnlyInRegion(const DeviceMesh& mesh,
                                             const std::string& contactName,
                                             const std::string& regionName)
{
    const auto nodeCells = buildNodeCellMap(mesh);
    std::vector<Index> nodes;

    for (const Contact& contact : mesh.contacts()) {
        if (contact.name != contactName)
            continue;

        for (Index nodeId : contact.node_ids) {
            bool hasCell = false;
            bool onlyRegion = true;
            for (Index cellId : nodeCells.at(nodeId)) {
                hasCell = true;
                const Region& region = mesh.getRegion(mesh.getCell(cellId).region_id);
                if (region.name != regionName) {
                    onlyRegion = false;
                    break;
                }
            }
            if (hasCell && onlyRegion)
                nodes.push_back(nodeId);
        }
        break;
    }

    return nodes;
}

Index firstSemiconductorOxideInterfaceEdge(const DeviceMesh& mesh,
                                           const MaterialDatabase& matdb)
{
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        bool hasTransportCell = false;
        bool hasInsulatingCell = false;
        for (Index cellId : edgeCells.at(edgeId)) {
            const Region& region = mesh.getRegion(mesh.getCell(cellId).region_id);
            const Material material = matdb.getMaterial(region.material);
            hasTransportCell =
                hasTransportCell || material.mun > 0.0 || material.mup > 0.0;
            hasInsulatingCell =
                hasInsulatingCell || (material.mun <= 0.0 && material.mup <= 0.0);
        }
        if (hasTransportCell && hasInsulatingCell)
            return edgeId;
    }

    FAIL("expected the MOS mesh to contain a semiconductor/oxide interface edge");
    return 0;
}

std::vector<std::filesystem::path> mosExampleDirs()
{
    const std::filesystem::path examplesRoot =
        std::filesystem::path(VELA_SOURCE_DIR) / "examples";
    return {
        examplesRoot / "nmos2d_mos_dd",
        examplesRoot / "pmos2d_mos_dd",
    };
}

DeviceMesh makeHorizontalInterfaceMesh()
{
    DeviceMesh mesh;
    const Real L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;  mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;  mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5 * L; n2.y = -L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.5 * L; n3.y = L;    mesh.addNode(n3);

    Cell si; si.id = 0; si.type = CellType::Tri3; si.region_id = 0;
    si.node_ids = {0, 1, 2}; mesh.addCell(si);
    Cell oxide; oxide.id = 1; oxide.type = CellType::Tri3; oxide.region_id = 1;
    oxide.node_ids = {0, 3, 1}; mesh.addCell(oxide);

    Region channel; channel.id = 0; channel.name = "channel"; channel.material = "Si";
    channel.cell_ids = {0}; mesh.addRegion(channel);
    Region gateOxide; gateOxide.id = 1; gateOxide.name = "gate_oxide"; gateOxide.material = "SiO2";
    gateOxide.cell_ids = {1}; mesh.addRegion(gateOxide);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeOxideFirstInterfaceMesh()
{
    DeviceMesh mesh;
    const Real L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;  mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;  mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5 * L; n2.y = L;    mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.5 * L; n3.y = -L;   mesh.addNode(n3);

    Cell oxide; oxide.id = 0; oxide.type = CellType::Tri3; oxide.region_id = 0;
    oxide.node_ids = {0, 1, 2}; mesh.addCell(oxide);
    Cell silicon; silicon.id = 1; silicon.type = CellType::Tri3; silicon.region_id = 1;
    silicon.node_ids = {0, 3, 1}; mesh.addCell(silicon);

    Region gateOxide; gateOxide.id = 0; gateOxide.name = "gate_oxide";
    gateOxide.material = "SiO2"; gateOxide.cell_ids = {0}; mesh.addRegion(gateOxide);
    Region channel; channel.id = 1; channel.name = "channel";
    channel.material = "Si"; channel.cell_ids = {1}; mesh.addRegion(channel);

    mesh.buildEdges();
    return mesh;
}

Index findEdgeByNodes(const DeviceMesh& mesh, Index a, Index b)
{
    if (b < a)
        std::swap(a, b);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        Index e0 = edge.n0;
        Index e1 = edge.n1;
        if (e1 < e0)
            std::swap(e0, e1);
        if (e0 == a && e1 == b)
            return edgeId;
    }
    FAIL("expected edge in test mesh");
    return 0;
}

} // namespace

TEST_CASE("interface intrinsic density prefers silicon when oxide cells are ordered first",
          "[mos_mixed][material]")
{
    const DeviceMesh mesh = makeOxideFirstInterfaceMesh();
    const MaterialDatabase matdb;
    const auto ni = detail::buildNodeNi(mesh, matdb, constants::T0);

    REQUIRE(ni.at(0) == matdb.getMaterial("Si").ni);
    REQUIRE(ni.at(1) == matdb.getMaterial("Si").ni);
    REQUIRE(ni.at(2) == 0.0);
    REQUIRE(ni.at(3) == matdb.getMaterial("Si").ni);
}

TEST_CASE("triangle GSS avalanche excludes oxide cells at a shared interface",
          "[mos_mixed][impact][triangle_gss]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    const auto materials = detail::buildCellMaterials(
        mesh, matdb, constants::T0);

    MobilityModelConfig mobilityConfig = mobilityModelConfig("masetti_field");
    mobilityConfig.highFieldDrivingForce = "quasi_fermi_gradient";
    const auto mobility = makeMobilityModel(mobilityConfig);
    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "selberherr";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "cell_reconstructed";
    impactConfig.sourceMappingMode = "triangle_gss_gradqf_truncated";
    impactConfig.sourceGeometryScale = 1.0;
    impactConfig.electronA = 1.0;
    impactConfig.electronB = 1.0e-30;
    impactConfig.holeA = 1.0;
    impactConfig.holeB = 1.0e-30;
    const auto impact = makeImpactIonizationModel(impactConfig);

    VectorXd psi(4);
    psi << 0.0, 0.1, -0.5, 0.8;
    VectorXd phin(4);
    phin << 0.0, -0.2, -0.6, -0.4;
    VectorXd phip(4);
    phip << -0.1, -0.3, -0.8, -0.5;
    const VectorXd n = VectorXd::Constant(4, 1.0e16);
    const VectorXd p = VectorXd::Constant(4, 1.0e16);

    const auto siliconRecords = detail::triangleGssAvalancheSourceRecordsForCell(
        impactConfig, *impact, mobilityConfig, *mobility, cellEdges.at(0),
        mesh, doping, materials, 0, psi, phin, phip, n, p,
        constants::Vt_300);
    const auto oxideRecords = detail::triangleGssAvalancheSourceRecordsForCell(
        impactConfig, *impact, mobilityConfig, *mobility, cellEdges.at(1),
        mesh, doping, materials, 1, psi, phin, phip, n, p,
        constants::Vt_300);

    REQUIRE(siliconRecords.size() == 3);
    REQUIRE(oxideRecords.empty());
}


TEST_CASE("surface mobility uses the reconstructed normal interface field",
          "[mos_mixed][mobility][surface]")
{
    DeviceMesh mesh = makeHorizontalInterfaceMesh();
    MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 0.0, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const Real temperature_K = constants::Vt_300 * constants::q / constants::kb;
    const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);

    MobilityModelConfig config = mobilityModelConfig("caughey_thomas_surface");
    config.surface.thetaElectron = 2.0e-6;
    config.surface.surfaceRegion = "channel";
    config.surface.surfaceInterface = {"channel", "gate_oxide"};
    const auto mobility = makeMobilityModel(config);

    VectorXd psi(4);
    psi << 0.0, 0.0, 0.0, 1.0;
    const Index interfaceEdge = findEdgeByNodes(mesh, 0, 1);

    const Real tangentialOnly = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config);
    const Real normalLimited = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi);

    REQUIRE(tangentialOnly > 0.0);
    REQUIRE(normalLimited > 0.0);
    REQUIRE(normalLimited < tangentialOnly);
}

TEST_CASE("surface normal field respects TCAD coordinate and field units",
          "[mos_mixed][mobility][surface][unit_scaling]")
{
    DeviceMesh mesh = makeHorizontalInterfaceMesh();
    MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 0.0, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(
        mesh, matdb, constants::T0);

    MobilityModelConfig config = mobilityModelConfigFromJson(
        nlohmann::json{
            {"model", "masetti_field_surface"},
            {"surface", {
                {"theta_electron_m_per_V", 1.0e-6},
                {"surface_region", "channel"},
                {"surface_interface", {"channel", "gate_oxide"}},
            }},
        },
        UnitScalingConfig{UnitScalingMode::UnitScaling});
    const auto mobility = makeMobilityModel(config);

    VectorXd psi(4);
    psi << 0.0, 0.0, 0.0, 1.0;
    const Index interfaceEdge = findEdgeByNodes(mesh, 0, 1);
    const Real withoutNormalField = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config);
    const Real withNormalField = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi);

    REQUIRE(config.surface.coordinateFieldFactor == 1.0e4);
    REQUIRE(withNormalField < withoutNormalField * 0.01);
}

TEST_CASE("Enhanced Lombardi geometry reaches cells beyond the interface row",
          "[mos_mixed][mobility][surface][lombardi][distance]")
{
    DeviceMesh mesh = makeHorizontalInterfaceMesh();
    MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto materials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    MobilityModelConfig config = mobilityModelConfig("masetti_lombardi");
    config.surface.surfaceRegion = "channel";
    config.surface.surfaceInterface = {"channel", "gate_oxide"};
    VectorXd psi(4);
    psi << 0.0, 0.0, -1.0, 1.0;

    detail::updateSurfaceMobilityCellGeometry(
        config, mesh, edgeCells, psi, 1.0);
    REQUIRE(config.surface.cellNormalFields.size() == mesh.numCells());
    REQUIRE(config.surface.cellDistances.size() == mesh.numCells());
    REQUIRE(std::isfinite(config.surface.cellNormalFields.at(0)));
    REQUIRE(std::isfinite(config.surface.cellDistances.at(0)));
    REQUIRE(config.surface.cellDistances.at(0) > 0.0);

    const auto mobility = makeMobilityModel(config);
    const Index interfaceEdge = findEdgeByNodes(mesh, 0, 1);
    const Real limited = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, materials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi);
    MobilityModelConfig bulkConfig = mobilityModelConfig("masetti");
    const auto bulk = makeMobilityModel(bulkConfig);
    const Real baseline = detail::edgeMobility(
        edgeCells, mesh, doping, *bulk, materials, interfaceEdge,
        CarrierType::Electron, 0.0, &bulkConfig, &psi);
    REQUIRE(limited < baseline);
}

TEST_CASE("mixed Si/SiO2 MOS edge mobility preserves semiconductor interface transport",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            const auto edgeCells = detail::buildEdgeCellMap(mesh);
            const Real temperature_K = constants::Vt_300 * constants::q / constants::kb;
            const auto cellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);
            const ConstantMobility mobility;
            const Index interfaceEdge = firstSemiconductorOxideInterfaceEdge(mesh, matdb);

            const Real mun = detail::edgeMobility(edgeCells,
                                                  mesh,
                                                  doping,
                                                  mobility,
                                                  cellMaterials,
                                                  interfaceEdge,
                                                  CarrierType::Electron,
                                                  0.0);
            const Real mup = detail::edgeMobility(edgeCells,
                                                  mesh,
                                                  doping,
                                                  mobility,
                                                  cellMaterials,
                                                  interfaceEdge,
                                                  CarrierType::Hole,
                                                  0.0);

            REQUIRE(mun > 0.0);
            REQUIRE(mup > 0.0);
        }
    }
}

TEST_CASE("mixed Si/SiO2 MOS DD scalar assembly keeps oxide carrier rows finite",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            DDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

            const int nNodes = static_cast<int>(mesh.numNodes());
            VectorXd psi = VectorXd::Zero(nNodes);
            VectorXd n = VectorXd::Constant(nNodes, 1.0e16);
            VectorXd p = VectorXd::Constant(nNodes, 1.0e16);

            assembler.assembleElectronContinuity(psi, n, p);
            REQUIRE(vectorIsFinite(assembler.rhs()));
            REQUIRE(assembler.matrix().rows() == nNodes);
            REQUIRE(assembler.matrix().cols() == nNodes);

            const std::vector<Index> oxideGateOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideGateOnlyNodes.empty());
            for (Index node : oxideGateOnlyNodes) {
                INFO("electron oxide node " << node);
                const int row = static_cast<int>(node);
                REQUIRE(assembler.matrix().coeff(row, row) == 1.0);
                REQUIRE(assembler.rhs()(row) == 0.0);
            }

            assembler.assembleHoleContinuity(psi, n, p);
            REQUIRE(vectorIsFinite(assembler.rhs()));
            for (Index node : oxideGateOnlyNodes) {
                INFO("hole oxide node " << node);
                const int row = static_cast<int>(node);
                REQUIRE(assembler.matrix().coeff(row, row) == 1.0);
                REQUIRE(assembler.rhs()(row) == 0.0);
            }
        }
    }
}

TEST_CASE("mixed Si/SiO2 MOS coupled DD residual and Jacobian are finite",
          "[mos_mixed][dd]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);

            const int nNodes = static_cast<int>(mesh.numNodes());
            CoupledDDState state;
            state.psi = VectorXd::Zero(nNodes);
            state.phin = VectorXd::Constant(nNodes, 0.25);
            state.phip = VectorXd::Constant(nNodes, -0.25);
            const VectorXd x = assembler.pack(state);

            // Pure oxide nodes legitimately have ni = n = p = 0.  They are
            // pinned algebraic rows, not invalid semiconductor carrier states.
            REQUIRE(assembler.hasPositiveFiniteCarriers(x));

            CoupledDDBoundaryConditions bcs;
            const VectorXd residual = assembler.residual(x, bcs);
            REQUIRE(residual.size() == 3 * nNodes);
            REQUIRE(vectorIsFinite(residual));

            const SparseMatrixd jacobian = assembler.assembleJacobian(x, bcs);
            REQUIRE(jacobian.rows() == 3 * nNodes);
            REQUIRE(jacobian.cols() == 3 * nNodes);
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator it(jacobian, outer); it; ++it)
                    REQUIRE(std::isfinite(it.value()));
            }

            const std::vector<Index> oxideGateOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideGateOnlyNodes.empty());
            for (Index node : oxideGateOnlyNodes) {
                INFO("coupled oxide node " << node);
                const int electronRow = nNodes + static_cast<int>(node);
                const int holeRow = 2 * nNodes + static_cast<int>(node);
                REQUIRE(residual(electronRow) == state.phin(static_cast<int>(node)));
                REQUIRE(residual(holeRow) == state.phip(static_cast<int>(node)));
                REQUIRE(jacobian.coeff(electronRow, electronRow) == 1.0);
                REQUIRE(jacobian.coeff(holeRow, holeRow) == 1.0);
            }
        }
    }
}

TEST_CASE("mixed Si/SiO2 MOS carrier gauges are quasi-Fermi reference invariant",
          "[mos_mixed][dd][qf-reference]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg =
                readJson(exampleDir / "simulation_iv.json");
            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            CoupledDDAssembler absolute(
                mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);
            CoupledDDAssembler referenced(
                mesh, matdb, doping, constants::Vt_300, 1.0e-6, 1.0e-6);
            referenced.setQuasiFermiReferences(1.1, -0.7);

            const int nNodes = static_cast<int>(mesh.numNodes());
            CoupledDDState state;
            state.psi = VectorXd::LinSpaced(nNodes, -0.15, 0.35);
            state.phin = VectorXd::LinSpaced(
                nNodes, 1.1 - 2.0e-4, 1.1 + 3.0e-4);
            state.phip = VectorXd::LinSpaced(
                nNodes, -0.7 + 4.0e-4, -0.7 - 1.0e-4);

            const VectorXd absoluteX = absolute.pack(state);
            const VectorXd referencedX = referenced.pack(state);
            const CoupledDDBoundaryConditions bcs;
            const VectorXd absoluteResidual = absolute.residual(absoluteX, bcs);
            const VectorXd referencedResidual =
                referenced.residual(referencedX, bcs);
            const Real residualScale =
                std::max<Real>(1.0, absoluteResidual.norm());
            REQUIRE((referencedResidual - absoluteResidual).norm()
                    <= 5.0e-12 * residualScale);

            const Eigen::MatrixXd absoluteJacobian(
                absolute.assembleJacobian(absoluteX, bcs));
            const Eigen::MatrixXd referencedJacobian(
                referenced.assembleJacobian(referencedX, bcs));
            const Real jacobianScale =
                std::max<Real>(1.0, absoluteJacobian.norm());
            Eigen::Index maxDifferenceRow = 0;
            Eigen::Index maxDifferenceColumn = 0;
            const Real maxJacobianDifference =
                (referencedJacobian - absoluteJacobian).cwiseAbs().maxCoeff(
                    &maxDifferenceRow, &maxDifferenceColumn);
            CAPTURE(maxDifferenceRow, maxDifferenceColumn,
                    maxJacobianDifference, jacobianScale);
            REQUIRE((referencedJacobian - absoluteJacobian).norm()
                    <= 5.0e-12 * jacobianScale);

            const auto absoluteTerms =
                absolute.carrierContinuityTermDiagnostics(absoluteX, bcs);
            const auto referencedTerms =
                referenced.carrierContinuityTermDiagnostics(referencedX, bcs);
            REQUIRE(referencedTerms.size() == absoluteTerms.size());
            for (std::size_t i = 0; i < absoluteTerms.size(); ++i) {
                CAPTURE(i);
                REQUIRE(referencedTerms[i].electronResidual ==
                        Catch::Approx(absoluteTerms[i].electronResidual)
                            .epsilon(5.0e-12).margin(1.0e-12));
                REQUIRE(referencedTerms[i].holeResidual ==
                        Catch::Approx(absoluteTerms[i].holeResidual)
                            .epsilon(5.0e-12).margin(1.0e-12));
            }

            const std::vector<Index> oxideGateOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideGateOnlyNodes.empty());
            for (Index node : oxideGateOnlyNodes) {
                const int ii = static_cast<int>(node);
                REQUIRE(referencedTerms[node].electronGauge ==
                        Catch::Approx(state.phin(ii)).margin(1.0e-15));
                REQUIRE(referencedTerms[node].holeGauge ==
                        Catch::Approx(state.phip(ii)).margin(1.0e-15));
            }
        }
    }
}

TEST_CASE("local-AD avalanche keeps pure oxide carrier rows as unit gauges",
          "[mos_mixed][dd][impact][local_ad]")
{
    for (const std::filesystem::path& exampleDir : mosExampleDirs()) {
        DYNAMIC_SECTION(exampleDir.filename().string()) {
            const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");
            JsonMeshReader reader;
            DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
            MaterialDatabase matdb;
            DopingModel doping = dopingFromDeck(mesh, cfg);

            MobilityModelConfig mobility = mobilityModelConfig("masetti_field");
            mobility.highFieldDrivingForce = "quasi_fermi_gradient";
            RecombinationModelConfig recombination;
            recombination.mechanisms = {"none"};
            ImpactIonizationModelConfig impact;
            impact.model = "selberherr";
            impact.drivingForce = "quasi_fermi_gradient";
            impact.generation = "current_density";
            impact.currentApproximation = "cell_reconstructed";
            impact.currentMagnitudeMode = "edge_scalar_abs";
            impact.cellReconstructedMidpointDensity = "gss_logistic";
            impact.quasiFermiGradientDiscretization = "cell_gradient";
            impact.sourceVolumePolicy = "genius_truncated";
            impact.sourceVolumeFactor = 0.0;
            impact.sourceGeometryScale = 1.0;
            impact.edgeSourcePartition = "symmetric";
            impact.sourceMappingMode = "triangle_gss_gradqf_truncated";
            impact.sourceJacobianMode = "local_ad";
            impact.electronA = 1.0;
            impact.electronB = 1.0e-30;
            impact.holeA = 1.0;
            impact.holeB = 1.0e-30;

            CoupledDDAssembler assembler(
                mesh, matdb, doping, constants::Vt_300, mobility,
                recombination, {}, impact);
            const int nNodes = static_cast<int>(mesh.numNodes());
            CoupledDDState state;
            state.psi = VectorXd::LinSpaced(nNodes, 0.2, -2.0);
            state.phin = VectorXd::LinSpaced(nNodes, 0.0, -1.2);
            state.phip = VectorXd::LinSpaced(nNodes, -0.1, -0.8);
            const VectorXd x = assembler.pack(state);
            const CoupledDDBoundaryConditions bcs;
            const VectorXd residual = assembler.residual(x, bcs);
            const SparseMatrixd jacobian = assembler.assembleJacobian(x, bcs);

            const std::vector<Index> oxideOnlyNodes =
                contactNodesOnlyInRegion(mesh, "gate", "gate_oxide");
            REQUIRE_FALSE(oxideOnlyNodes.empty());
            for (Index node : oxideOnlyNodes) {
                for (const int row : {
                         nNodes + static_cast<int>(node),
                         2 * nNodes + static_cast<int>(node)}) {
                    INFO("pure oxide carrier row " << row);
                    int numericEntries = 0;
                    for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                        for (SparseMatrixd::InnerIterator entry(jacobian, outer);
                             entry; ++entry) {
                            if (entry.row() != row || entry.value() == 0.0)
                                continue;
                            ++numericEntries;
                            REQUIRE(entry.col() == row);
                            REQUIRE(entry.value() == 1.0);
                        }
                    }
                    REQUIRE(numericEntries == 1);
                    REQUIRE(residual(row) == x(row));
                }
            }
        }
    }
}
