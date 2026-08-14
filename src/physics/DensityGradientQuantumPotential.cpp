#include "vela/physics/DensityGradientQuantumPotential.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/solver/LinearSolver.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace vela {
namespace {

Real safeExponential(Real argument)
{
    constexpr Real limit = 700.0;
    return std::exp(std::clamp(argument, -limit, limit));
}

struct Eq231CellResidualRecord {
    Index cellId = 0;
    Index regionId = 0;
    bool isTransport = true;
    std::array<Index, 3> nodes{};
    std::array<Real, 3> lambda_V{};
    std::array<Real, 3> auxiliary{};
    std::array<Real, 3> stiffness{};
    std::array<Real, 3> gradientSquared{};
    std::array<Real, 3> reaction{};
    Real area_m2 = 0.0;
    Real gradientSquaredPerM2 = 0.0;
};

struct Eq231ResidualDiagnostic {
    std::vector<Eq231CellResidualRecord> cells;
    VectorXd stiffness;
    VectorXd gradientSquared;
    VectorXd reaction;
    VectorXd interfaceBoundary;
};

bool isSiliconDioxide(const std::string& material)
{
    std::string lower;
    lower.reserve(material.size());
    for (const unsigned char ch : material)
        lower.push_back(static_cast<char>(std::tolower(ch)));
    return lower == "sio2" || lower == "oxide" ||
           lower == "silicondioxide";
}

std::string csvField(const std::string& value)
{
    if (value.find_first_of(",\"\r\n") == std::string::npos)
        return value;
    std::string quoted = "\"";
    for (const char ch : value) {
        quoted += ch;
        if (ch == '"')
            quoted += '"';
    }
    quoted += '"';
    return quoted;
}

void writeEq231ResidualDiagnostic(
    const std::string& prefixValue,
    const DeviceMesh& mesh,
    const Eq231ResidualDiagnostic& diagnostic,
    const VectorXd& state,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V)
{
    namespace fs = std::filesystem;
    const fs::path prefix(prefixValue);
    if (!prefix.parent_path().empty())
        fs::create_directories(prefix.parent_path());
    const auto open = [&](const std::string& suffix) {
        std::ofstream stream(prefix.string() + suffix);
        if (!stream)
            throw std::runtime_error(
                "failed to create Eq. 231 residual diagnostic: " +
                prefix.string() + suffix);
        stream << std::setprecision(17);
        return stream;
    };

    using EdgeKey = std::pair<Index, Index>;
    std::map<EdgeKey, std::vector<Index>> edgeCells;
    for (const auto& record : diagnostic.cells) {
        for (int local = 0; local < 3; ++local) {
            Index a = record.nodes[local];
            Index b = record.nodes[(local + 1) % 3];
            if (a > b)
                std::swap(a, b);
            edgeCells[{a, b}].push_back(record.cellId);
        }
    }
    std::map<Index, std::set<std::pair<Index, Index>>> interfacePairsByCell;
    for (const auto& [edge, cells] : edgeCells) {
        (void)edge;
        if (cells.size() != 2)
            continue;
        const Index region0 = mesh.getCell(cells[0]).region_id;
        const Index region1 = mesh.getCell(cells[1]).region_id;
        if (region0 == region1)
            continue;
        const auto pair = std::minmax(region0, region1);
        interfacePairsByCell[cells[0]].insert(pair);
        interfacePairsByCell[cells[1]].insert(pair);
    }
    const auto pairText = [&](Index cellId) {
        std::ostringstream value;
        const auto found = interfacePairsByCell.find(cellId);
        if (found == interfacePairsByCell.end())
            return value.str();
        bool first = true;
        for (const auto& [region0, region1] : found->second) {
            if (!first)
                value << ';';
            first = false;
            value << mesh.getRegion(region0).name << '|'
                  << mesh.getRegion(region1).name;
        }
        return value.str();
    };

    auto cells = open("_cells.csv");
    cells << "cell_id,region_id,region_name,material,is_interface_cell,"
             "interface_pairs,is_transport,local_node,node_id,x_internal,y_internal,"
             "is_active,is_dirichlet,area_m2,lambda_V,w,grad_w_squared_per_m2,"
             "stiffness,gradient_squared,reaction,total\n";
    for (const auto& record : diagnostic.cells) {
        const Region& region = mesh.getRegion(record.regionId);
        const bool interfaceCell = interfacePairsByCell.contains(record.cellId);
        const std::string pairs = pairText(record.cellId);
        for (int local = 0; local < 3; ++local) {
            const Index node = record.nodes[local];
            const Node& point = mesh.getNode(node);
            const Real total = record.stiffness[local] +
                record.gradientSquared[local] + record.reaction[local];
            cells << record.cellId << ',' << record.regionId << ','
                  << csvField(region.name) << ',' << csvField(region.material) << ','
                  << (interfaceCell ? 1 : 0) << ',' << csvField(pairs) << ','
                  << (record.isTransport ? 1 : 0) << ','
                  << local << ',' << node << ',' << point.x << ',' << point.y << ','
                  << (activeNodes[node] ? 1 : 0) << ','
                  << (dirichletLambda_V.contains(node) ? 1 : 0) << ','
                  << record.area_m2 << ',' << record.lambda_V[local] << ','
                  << record.auxiliary[local] << ','
                  << record.gradientSquaredPerM2 << ','
                  << record.stiffness[local] << ','
                  << record.gradientSquared[local] << ','
                  << record.reaction[local] << ',' << total << '\n';
        }
    }

    auto nodes = open("_nodes.csv");
    nodes << "node_id,x_internal,y_internal,is_active,is_dirichlet,"
             "potential_like_V,output_lambda_V,stiffness,gradient_squared,"
             "reaction,raw_total,solver_total,interface_boundary\n";
    Index maxFreeNode = 0;
    Real maxFreeResidual = -1.0;
    for (int row = 0; row < state.size(); ++row) {
        const Index node = static_cast<Index>(row);
        const Node& point = mesh.getNode(node);
        const bool fixed = dirichletLambda_V.contains(node);
        const Real rawTotal = diagnostic.stiffness(row) +
            diagnostic.gradientSquared(row) + diagnostic.reaction(row) +
            diagnostic.interfaceBoundary(row);
        Real solverTotal = rawTotal;
        if (!activeNodes[node])
            solverTotal = state(row);
        if (fixed) {
            solverTotal = state(row) -
                (dirichletLambda_V.at(node) - nodeOutputBandDrive_V(row));
        }
        if (activeNodes[node] && !fixed && std::abs(rawTotal) > maxFreeResidual) {
            maxFreeResidual = std::abs(rawTotal);
            maxFreeNode = node;
        }
        nodes << node << ',' << point.x << ',' << point.y << ','
              << (activeNodes[node] ? 1 : 0) << ',' << (fixed ? 1 : 0) << ','
              << state(row) << ',' << state(row) + nodeOutputBandDrive_V(row) << ','
              << diagnostic.stiffness(row) << ','
              << diagnostic.gradientSquared(row) << ','
              << diagnostic.reaction(row) << ',' << rawTotal << ','
              << solverTotal << ',' << diagnostic.interfaceBoundary(row) << '\n';
    }

    struct RegionAggregate {
        std::set<Index> cells;
        std::set<Index> interfaceCells;
        Real stiffnessL1Free = 0.0;
        Real gradientL1Free = 0.0;
        Real reactionL1Free = 0.0;
        Real totalL1Free = 0.0;
        Real totalL1All = 0.0;
        Real interfaceTotalL1Free = 0.0;
        Real maxCellResidual = -1.0;
        Index maxCellId = 0;
    };
    std::map<Index, RegionAggregate> aggregates;
    Real globalStiffnessL1Free = 0.0;
    Real globalGradientL1Free = 0.0;
    Real globalReactionL1Free = 0.0;
    Real globalInterfaceBoundaryL1Free = 0.0;
    Real globalTotalL1Free = 0.0;
    for (const auto& record : diagnostic.cells) {
        auto& aggregate = aggregates[record.regionId];
        aggregate.cells.insert(record.cellId);
        const bool interfaceCell = interfacePairsByCell.contains(record.cellId);
        if (interfaceCell)
            aggregate.interfaceCells.insert(record.cellId);
        Real cellMaxFree = 0.0;
        for (int local = 0; local < 3; ++local) {
            const Index node = record.nodes[local];
            const Real total = record.stiffness[local] +
                record.gradientSquared[local] + record.reaction[local];
            aggregate.totalL1All += std::abs(total);
            if (!activeNodes[node] || dirichletLambda_V.contains(node))
                continue;
            aggregate.stiffnessL1Free += std::abs(record.stiffness[local]);
            aggregate.gradientL1Free += std::abs(record.gradientSquared[local]);
            aggregate.reactionL1Free += std::abs(record.reaction[local]);
            aggregate.totalL1Free += std::abs(total);
            globalStiffnessL1Free += std::abs(record.stiffness[local]);
            globalGradientL1Free += std::abs(record.gradientSquared[local]);
            globalReactionL1Free += std::abs(record.reaction[local]);
            globalTotalL1Free += std::abs(total);
            if (interfaceCell)
                aggregate.interfaceTotalL1Free += std::abs(total);
            cellMaxFree = std::max(cellMaxFree, std::abs(total));
        }
        if (cellMaxFree > aggregate.maxCellResidual) {
            aggregate.maxCellResidual = cellMaxFree;
            aggregate.maxCellId = record.cellId;
        }
    }
    for (int row = 0; row < diagnostic.interfaceBoundary.size(); ++row) {
        const Index node = static_cast<Index>(row);
        if (activeNodes[node] && !dirichletLambda_V.contains(node))
            globalInterfaceBoundaryL1Free +=
                std::abs(diagnostic.interfaceBoundary(row));
    }
    auto regions = open("_regions.csv");
    regions << "region_id,region_name,material,cell_count,interface_cell_count,"
               "stiffness_l1_free,gradient_squared_l1_free,reaction_l1_free,"
               "total_l1_free,total_l1_all,interface_total_l1_free,"
               "max_cell_residual_free,max_cell_id\n";
    for (const auto& [regionId, aggregate] : aggregates) {
        const Region& region = mesh.getRegion(regionId);
        regions << regionId << ',' << csvField(region.name) << ','
                << csvField(region.material) << ',' << aggregate.cells.size() << ','
                << aggregate.interfaceCells.size() << ','
                << aggregate.stiffnessL1Free << ','
                << aggregate.gradientL1Free << ','
                << aggregate.reactionL1Free << ','
                << aggregate.totalL1Free << ',' << aggregate.totalL1All << ','
                << aggregate.interfaceTotalL1Free << ','
                << aggregate.maxCellResidual << ',' << aggregate.maxCellId << '\n';
    }

    auto summary = open("_summary.txt");
    summary << "max_free_node=" << maxFreeNode << '\n'
            << "max_free_residual=" << maxFreeResidual << '\n'
            << "stiffness_l1_free=" << globalStiffnessL1Free << '\n'
            << "gradient_squared_l1_free=" << globalGradientL1Free << '\n'
            << "reaction_l1_free=" << globalReactionL1Free << '\n'
            << "interface_boundary_l1_free="
            << globalInterfaceBoundaryL1Free << '\n'
            << "cell_total_l1_free=" << globalTotalL1Free << '\n';
}

} // namespace

Real densityGradientStepFunction(Real x)
{
    if (!std::isfinite(x))
        throw std::invalid_argument("density-gradient step argument must be finite.");
    const Real ax = std::abs(x);
    if (ax < 1.0e-4) {
        // 1/2 + x/6 + x^2/24 + x^3/120 + x^4/720
        return 0.5 + x * (1.0 / 6.0 + x * (1.0 / 24.0 +
            x * (1.0 / 120.0 + x / 720.0)));
    }
    if (x > 700.0)
        return std::numeric_limits<Real>::max();
    return (std::expm1(x) - x) / (x * x);
}

Real densityGradientStepFunctionDerivative(Real x)
{
    if (!std::isfinite(x))
        throw std::invalid_argument("density-gradient step argument must be finite.");
    const Real ax = std::abs(x);
    if (ax < 1.0e-4) {
        // 1/6 + x/12 + x^2/40 + x^3/180 + x^4/1008
        return 1.0 / 6.0 + x * (1.0 / 12.0 + x * (1.0 / 40.0 +
            x * (1.0 / 180.0 + x / 1008.0)));
    }
    if (x > 700.0)
        return std::numeric_limits<Real>::max();
    return (std::exp(x) * (x - 2.0) + x + 2.0) / (x * x * x);
}

Real densityGradientCoefficientVm2(
    const DensityGradientQuantumPotentialConfig& config)
{
    return densityGradientCoefficientVm2(
        config.gamma, config.effectiveMassRatio);
}

Real densityGradientCoefficientVm2(Real gamma, Real effectiveMassRatio)
{
    if (!(gamma > 0.0) || !std::isfinite(gamma) ||
        !(effectiveMassRatio > 0.0) ||
        !std::isfinite(effectiveMassRatio)) {
        throw std::invalid_argument(
            "density-gradient gamma and effective mass must be finite and positive.");
    }
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    return gamma * hbar * hbar /
        (6.0 * effectiveMassRatio * constants::m0 * constants::q);
}

Real densityGradientOxidePenetrationDepthM(
    Real barrierEffectiveMassRatio,
    Real barrierHeight_V)
{
    if (!(barrierEffectiveMassRatio > 0.0) ||
        !std::isfinite(barrierEffectiveMassRatio) ||
        !(barrierHeight_V > 0.0) || !std::isfinite(barrierHeight_V)) {
        throw std::invalid_argument(
            "density-gradient oxide WKB mass and barrier must be positive.");
    }
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    return hbar / std::sqrt(
        2.0 * barrierEffectiveMassRatio * constants::m0 *
        constants::q * barrierHeight_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    const VectorXd coefficient = VectorXd::Constant(
        static_cast<int>(mesh.numNodes()), densityGradientCoefficientVm2(config));
    return solveElectronDensityGradientPotential(
        mesh, classicalDensity, coefficient, activeNodes, dirichletPotential_V, {},
        thermalVoltage_V, units, std::move(config), initialPotential_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    return solveElectronDensityGradientPotential(
        mesh, classicalDensity, coefficientVm2, activeNodes,
        dirichletPotential_V, {}, thermalVoltage_V, units,
        std::move(config), initialPotential_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    const std::vector<DensityGradientStepBoundary>& stepBoundaries,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    const int nodeCount = static_cast<int>(mesh.numNodes());
    if (classicalDensity.size() != nodeCount ||
        coefficientVm2.size() != nodeCount ||
        activeNodes.size() != static_cast<std::size_t>(nodeCount)) {
        throw std::invalid_argument(
            "density-gradient input sizes must match the mesh node count.");
    }
    if (!(thermalVoltage_V > 0.0) || !std::isfinite(thermalVoltage_V))
        throw std::invalid_argument("density-gradient thermal voltage must be positive.");
    if (config.maxIterations <= 0 || !(config.relativeTolerance >= 0.0) ||
        !(config.absoluteTolerance_V >= 0.0) || !(config.damping > 0.0) ||
        config.damping > 1.0 || !std::isfinite(config.damping) ||
        !(config.maxUpdate_V >= 0.0) || !std::isfinite(config.maxUpdate_V)) {
        throw std::invalid_argument("invalid density-gradient iteration configuration.");
    }

    VectorXd potential = initialPotential_V.size() == nodeCount
        ? initialPotential_V : VectorXd::Zero(nodeCount);
    for (int i = 0; i < nodeCount; ++i) {
        if (!std::isfinite(classicalDensity(i)) || classicalDensity(i) < 0.0)
            throw std::invalid_argument("density-gradient density must be finite and non-negative.");
        if (activeNodes[static_cast<std::size_t>(i)] &&
            (!(coefficientVm2(i) > 0.0) || !std::isfinite(coefficientVm2(i)))) {
            throw std::invalid_argument(
                "density-gradient material coefficient must be finite and positive.");
        }
        if (!std::isfinite(potential(i)))
            potential(i) = 0.0;
        if (!activeNodes[static_cast<std::size_t>(i)])
            potential(i) = 0.0;
    }
    for (const auto& [node, value] : dirichletPotential_V) {
        if (node >= mesh.numNodes() || !std::isfinite(value))
            throw std::invalid_argument("invalid density-gradient Dirichlet boundary.");
        potential(static_cast<int>(node)) = value;
    }
    for (const auto& boundary : stepBoundaries) {
        if (boundary.edgeId >= mesh.numEdges() ||
            !std::isfinite(boundary.barrierHeightN0_V) ||
            !std::isfinite(boundary.barrierHeightN1_V) ||
            !(boundary.barrierEffectiveMassRatio > 0.0) ||
            !std::isfinite(boundary.barrierEffectiveMassRatio) ||
            !(boundary.solvedGamma > 0.0) ||
            !std::isfinite(boundary.solvedGamma) ||
            !(boundary.theta > 0.0) || !std::isfinite(boundary.theta) ||
            !(boundary.alphaDeterminantCubeRoot > 0.0) ||
            !std::isfinite(boundary.alphaDeterminantCubeRoot) ||
            !std::isfinite(boundary.normalDrivingPotentialGradient_V_per_m)) {
            throw std::invalid_argument("invalid density-gradient step boundary.");
        }
        const Edge& edge = mesh.getEdge(boundary.edgeId);
        if (!activeNodes[edge.n0] || !activeNodes[edge.n1])
            throw std::invalid_argument(
                "density-gradient step boundary endpoints must be active.");
    }

    const Real lengthScale_m = units.lengthMPerInternal();
    const Real areaScale_m2 = lengthScale_m * lengthScale_m;
    const Real referenceDensity = std::max<Real>(classicalDensity.maxCoeff(), 1.0);
    const Real densityFloor = std::numeric_limits<Real>::min();
    LinearSolver solver;
    DensityGradientQuantumPotentialResult result;

    for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
        VectorXd rootDensity(nodeCount);
        VectorXd derivative(nodeCount);
        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                rootDensity(i) = 0.0;
                derivative(i) = 0.0;
                continue;
            }
            const Real normalized = std::max(classicalDensity(i), densityFloor) /
                referenceDensity;
            rootDensity(i) = std::sqrt(normalized) *
                safeExponential(-potential(i) / (2.0 * thermalVoltage_V));
            derivative(i) = -rootDensity(i) / (2.0 * thermalVoltage_V);
        }

        VectorXd residual = VectorXd::Zero(nodeCount);
        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(static_cast<std::size_t>(nodeCount) * 8);
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            if (!activeNodes[edge.n0] || !activeNodes[edge.n1] ||
                edge.length <= 1.0e-30 || edge.couple <= 0.0)
                continue;
            const Real geometryWeight = edge.couple / edge.length;
            const Real laplaceWeightI = coefficientVm2(i) * geometryWeight;
            const Real laplaceWeightJ = coefficientVm2(j) * geometryWeight;
            residual(i) += laplaceWeightI * (rootDensity(j) - rootDensity(i));
            residual(j) += laplaceWeightJ * (rootDensity(i) - rootDensity(j));
            triplets.emplace_back(i, i, -laplaceWeightI * derivative(i));
            triplets.emplace_back(i, j,  laplaceWeightI * derivative(j));
            triplets.emplace_back(j, i,  laplaceWeightJ * derivative(i));
            triplets.emplace_back(j, j, -laplaceWeightJ * derivative(j));
        }

        // Sentaurus Eq. 233: analytic 1-D step condition on an interface
        // between a solved semiconductor and an unsolved nonmetal region.
        // In the transformed u equation, C*dn(u) =
        // C*u*(dn(drive)-dn(Lambda))/(2*Vt).  Linear edge shape functions
        // give each endpoint half of the physical interface-segment length.
        const Real hbar = constants::h / (2.0 * std::acos(-1.0));
        const Real kT_J = thermalVoltage_V * constants::q;
        for (const auto& boundary : stepBoundaries) {
            const Edge& edge = mesh.getEdge(boundary.edgeId);
            const Real endpointLength_m = 0.5 * edge.length * lengthScale_m;
            const Real kappa = std::sqrt(
                24.0 * boundary.barrierEffectiveMassRatio * constants::m0 * kT_J /
                (hbar * hbar * boundary.solvedGamma *
                 boundary.alphaDeterminantCubeRoot));
            const int rows[2] = {
                static_cast<int>(edge.n0), static_cast<int>(edge.n1)};
            const Real barriers[2] = {
                boundary.barrierHeightN0_V, boundary.barrierHeightN1_V};
            for (int endpoint = 0; endpoint < 2; ++endpoint) {
                const int row = rows[endpoint];
                const Real lambda = potential(row);
                const Real barrier = barriers[endpoint];
                const Real argument =
                    2.0 * boundary.theta * (lambda - barrier) / thermalVoltage_V;
                const Real step = densityGradientStepFunction(argument);
                const Real stepDerivative =
                    densityGradientStepFunctionDerivative(argument);
                const Real normalLambdaGradient =
                    kappa * (barrier - lambda) * step;
                const Real normalLambdaGradientDerivative = -kappa * (
                    step + (lambda - barrier) *
                        (2.0 * boundary.theta / thermalVoltage_V) *
                        stepDerivative);
                const Real scale = coefficientVm2(row) * endpointLength_m /
                    (2.0 * thermalVoltage_V);
                const Real transformedNormalGradient =
                    boundary.normalDrivingPotentialGradient_V_per_m -
                    normalLambdaGradient;
                residual(row) += scale * rootDensity(row) * transformedNormalGradient;
                triplets.emplace_back(
                    row, row,
                    scale * (derivative(row) * transformedNormalGradient -
                        rootDensity(row) * normalLambdaGradientDerivative));
            }
        }

        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                residual(i) = potential(i);
                triplets.emplace_back(i, i, 1.0);
                continue;
            }
            const Real volume_m2 = mesh.getNode(static_cast<Index>(i)).volume * areaScale_m2;
            residual(i) += potential(i) * rootDensity(i) * volume_m2;
            triplets.emplace_back(
                i, i, (rootDensity(i) + potential(i) * derivative(i)) * volume_m2);
        }

        for (const auto& [node, value] : dirichletPotential_V) {
            const int row = static_cast<int>(node);
            residual(row) = potential(row) - value;
            triplets.emplace_back(row, row, 1.0);
        }

        SparseMatrixd jacobian(nodeCount, nodeCount);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
        // Replace constrained rows after summation, including edge entries.
        for (const auto& [node, value] : dirichletPotential_V) {
            (void)value;
            const int row = static_cast<int>(node);
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        for (int row = 0; row < nodeCount; ++row) {
            if (activeNodes[static_cast<std::size_t>(row)])
                continue;
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        jacobian.prune(0.0);

        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        VectorXd update = solver.solve(jacobian, -residual);
        Real updateNorm = update.lpNorm<Eigen::Infinity>();
        Eigen::Index maxUpdateIndex = 0;
        update.cwiseAbs().maxCoeff(&maxUpdateIndex);
        if (config.maxUpdate_V > 0.0 && updateNorm > config.maxUpdate_V) {
            update *= config.maxUpdate_V / updateNorm;
            updateNorm = config.maxUpdate_V;
        }
        potential += config.damping * update;
        for (const auto& [node, value] : dirichletPotential_V)
            potential(static_cast<int>(node)) = value;
        result.iterations = iteration;
        result.lastUpdateInfinityNorm_V = config.damping * updateNorm;
        result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
        result.maxUpdateNodeValue_V = config.damping * update(maxUpdateIndex);
        const Real scale = std::max<Real>(potential.lpNorm<Eigen::Infinity>(), 1.0);
        if (updateNorm <= config.absoluteTolerance_V + config.relativeTolerance * scale) {
            result.converged = true;
            break;
        }
    }
    result.potentialInfinityNorm_V = potential.lpNorm<Eigen::Infinity>();
    result.potential_V = std::move(potential);
    return result;
}

DensityGradientQuantumPotentialResult
solveElectronDensityGradientPotentialLikeGlobal(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialLambda_V,
    const VectorXd& initialPotentialLike_V)
{
    const int nodeCount = static_cast<int>(mesh.numNodes());
    if (nodeOutputBandDrive_V.size() != nodeCount ||
        activeNodes.size() != static_cast<std::size_t>(nodeCount) ||
        !(thermalVoltage_V > 0.0) || !std::isfinite(thermalVoltage_V)) {
        throw std::invalid_argument("invalid global density-gradient input.");
    }
    if (config.globalDiscretization != "p1_direct" &&
        config.globalDiscretization != "exponential_fitted" &&
        config.globalDiscretization != "cvfem_full") {
        throw std::invalid_argument(
            "unsupported global density-gradient discretization.");
    }
    if (config.oxideBoundary != "none" &&
        config.oxideBoundary != "devsim_wkb") {
        throw std::invalid_argument(
            "unsupported global density-gradient oxide boundary.");
    }
    if (config.oxideBoundary == "devsim_wkb" &&
        (!(config.oxideQuantumMassRatio > 0.0) ||
         !(config.oxideBarrierMassRatio > 0.0) ||
         !(config.oxideBarrierHeight_V > 0.0) ||
         !std::isfinite(config.oxideQuantumMassRatio) ||
         !std::isfinite(config.oxideBarrierMassRatio) ||
         !std::isfinite(config.oxideBarrierHeight_V))) {
        throw std::invalid_argument(
            "invalid global density-gradient oxide boundary parameters.");
    }
    VectorXd continuousPotential = -nodeOutputBandDrive_V;
    if (initialPotentialLike_V.size() == nodeCount)
        continuousPotential = initialPotentialLike_V;
    else if (initialLambda_V.size() == nodeCount)
        continuousPotential = initialLambda_V - nodeOutputBandDrive_V;
    for (int i = 0; i < nodeCount; ++i) {
        if (!std::isfinite(continuousPotential(i)))
            continuousPotential(i) = -nodeOutputBandDrive_V(i);
        if (!activeNodes[static_cast<std::size_t>(i)])
            continuousPotential(i) = 0.0;
    }
    for (const auto& [node, lambda] : dirichletLambda_V) {
        if (node >= mesh.numNodes() || !std::isfinite(lambda))
            throw std::invalid_argument("invalid global density-gradient Dirichlet boundary.");
        continuousPotential(static_cast<int>(node)) =
            lambda - nodeOutputBandDrive_V(static_cast<int>(node));
    }
    for (const auto& data : cellMaterials) {
        if (data.cellId >= mesh.numCells() || !(data.coefficientVm2 > 0.0) ||
            !std::isfinite(data.coefficientVm2) ||
            mesh.getCell(data.cellId).node_ids.size() != 3) {
            throw std::invalid_argument("invalid global density-gradient cell material.");
        }
        for (int local = 0; local < 3; ++local) {
            if (!std::isfinite(data.materialBandDrive_V[local]) ||
                !std::isfinite(data.dynamicDrivingPotential_V[local]) ||
                !std::isfinite(data.initialLambda_V[local]))
                throw std::invalid_argument("nonfinite global density-gradient cell drive.");
        }
    }

    struct OxideInterfaceSegment {
        Index edgeId = 0;
        Index oxideCellId = 0;
    };
    std::vector<OxideInterfaceSegment> oxideInterfaceSegments;
    const bool useDevsimOxideBoundary = config.oxideBoundary == "devsim_wkb";
    const Real devsimOxideCoefficientVm2 = useDevsimOxideBoundary
        ? densityGradientCoefficientVm2(1.0, config.oxideQuantumMassRatio)
        : 0.0;
    const Real oxidePenetrationDepth_m = useDevsimOxideBoundary
        ? densityGradientOxidePenetrationDepthM(
            config.oxideBarrierMassRatio, config.oxideBarrierHeight_V)
        : 0.0;
    if (useDevsimOxideBoundary) {
        using EdgeKey = std::pair<Index, Index>;
        struct EdgeSide {
            Index cellId = 0;
            bool isTransport = false;
            bool isOxide = false;
        };
        std::map<EdgeKey, std::vector<EdgeSide>> sides;
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            const Region& region = mesh.getRegion(cell.region_id);
            for (int local = 0; local < 3; ++local) {
                Index a = cell.node_ids[local];
                Index b = cell.node_ids[(local + 1) % 3];
                if (a > b)
                    std::swap(a, b);
                sides[{a, b}].push_back({
                    data.cellId, data.isTransport,
                    isSiliconDioxide(region.material)});
            }
        }
        std::map<EdgeKey, Index> edgeIds;
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            edgeIds[std::minmax(edge.n0, edge.n1)] = edgeId;
        }
        for (const auto& [key, adjacent] : sides) {
            if (adjacent.size() != 2)
                continue;
            const EdgeSide* oxide = nullptr;
            const EdgeSide* transport = nullptr;
            for (const auto& side : adjacent) {
                if (!side.isTransport && side.isOxide)
                    oxide = &side;
                if (side.isTransport)
                    transport = &side;
            }
            if (oxide != nullptr && transport != nullptr)
                oxideInterfaceSegments.push_back({edgeIds.at(key), oxide->cellId});
        }
    }

    const Real lengthScale_m = units.lengthMPerInternal();
    const Real areaScale_m2 = lengthScale_m * lengthScale_m;
    LinearSolver solver;
    DensityGradientQuantumPotentialResult result;
    Real initialResidualNorm = -1.0;
    auto assemble = [&](const VectorXd& state,
                        VectorXd& residual,
                        std::vector<Eigen::Triplet<Real>>* triplets,
                        Eq231ResidualDiagnostic* diagnostic) {
        residual = VectorXd::Zero(nodeCount);
        if (triplets != nullptr) {
            triplets->clear();
                triplets->reserve(cellMaterials.size() * 18 + nodeCount);
        }
        if (diagnostic != nullptr) {
            diagnostic->cells.clear();
            diagnostic->cells.reserve(cellMaterials.size());
            diagnostic->stiffness = VectorXd::Zero(nodeCount);
            diagnostic->gradientSquared = VectorXd::Zero(nodeCount);
            diagnostic->reaction = VectorXd::Zero(nodeCount);
            diagnostic->interfaceBoundary = VectorXd::Zero(nodeCount);
        }
        auto addJacobian = [&](int row, int col, Real value) {
            if (triplets != nullptr)
                triplets->emplace_back(row, col, value);
        };
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            const Index nodes[3] = {
                cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
            const Node& p0 = mesh.getNode(nodes[0]);
            const Node& p1 = mesh.getNode(nodes[1]);
            const Node& p2 = mesh.getNode(nodes[2]);
            const Real twiceAreaInternal = std::abs(
                (p1.x - p0.x) * (p2.y - p0.y) -
                (p1.y - p0.y) * (p2.x - p0.x));
            if (!(twiceAreaInternal > 0.0))
                continue;
            const Real cellArea_m2 = 0.5 * twiceAreaInternal * areaScale_m2;
            const Region& cellRegion = mesh.getRegion(cell.region_id);
            const Real cellCoefficientVm2 =
                useDevsimOxideBoundary &&
                    isSiliconDioxide(cellRegion.material)
                ? devsimOxideCoefficientVm2
                : data.coefficientVm2;
            Real lambda[3]{};
            Real w[3]{};
            for (int local = 0; local < 3; ++local) {
                const int node = static_cast<int>(nodes[local]);
                lambda[local] = state(node) + data.materialBandDrive_V[local];
                w[local] = (data.dynamicDrivingPotential_V[local] -
                    state(node)) / thermalVoltage_V;
            }
            const Real x[3] = {p0.x * lengthScale_m,
                               p1.x * lengthScale_m,
                               p2.x * lengthScale_m};
            const Real y[3] = {p0.y * lengthScale_m,
                               p1.y * lengthScale_m,
                               p2.y * lengthScale_m};
            const Real b[3] = {y[1] - y[2], y[2] - y[0], y[0] - y[1]};
            const Real c[3] = {x[2] - x[1], x[0] - x[2], x[1] - x[0]};
            const Real fourArea2 = 4.0 * cellArea_m2 * cellArea_m2;
            Real cvfemFlux[3][3]{};
            if (config.globalDiscretization == "cvfem_full") {
                const Real centroidX = (x[0] + x[1] + x[2]) / 3.0;
                const Real centroidY = (y[0] + y[1] + y[2]) / 3.0;
                constexpr int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& pair : pairs) {
                    const int a = pair[0];
                    const int neighbour = pair[1];
                    const Real midpointX = 0.5 * (x[a] + x[neighbour]);
                    const Real midpointY = 0.5 * (y[a] + y[neighbour]);
                    const Real tangentX = centroidX - midpointX;
                    const Real tangentY = centroidY - midpointY;
                    Real normalLengthX = tangentY;
                    Real normalLengthY = -tangentX;
                    const Real towardNeighbourX = x[neighbour] - x[a];
                    const Real towardNeighbourY = y[neighbour] - y[a];
                    if (normalLengthX * towardNeighbourX +
                            normalLengthY * towardNeighbourY < 0.0) {
                        normalLengthX = -normalLengthX;
                        normalLengthY = -normalLengthY;
                    }
                    for (int d = 0; d < 3; ++d) {
                        const Real coefficient =
                            (b[d] * normalLengthX + c[d] * normalLengthY) /
                            (2.0 * cellArea_m2);
                        cvfemFlux[a][d] += coefficient;
                        cvfemFlux[neighbour][d] -= coefficient;
                    }
                }
            }
            Real gradientSquared = 0.0;
            for (int d = 0; d < 3; ++d) {
                for (int e = 0; e < 3; ++e) {
                    gradientSquared += w[d] * w[e] *
                        (b[d] * b[e] + c[d] * c[e]) / fourArea2;
                }
            }
            Eq231CellResidualRecord record;
            if (diagnostic != nullptr) {
                record.cellId = data.cellId;
                record.regionId = cell.region_id;
                record.isTransport = data.isTransport;
                record.nodes = {nodes[0], nodes[1], nodes[2]};
                record.area_m2 = cellArea_m2;
                record.gradientSquaredPerM2 = gradientSquared;
                for (int local = 0; local < 3; ++local) {
                    record.lambda_V[local] = lambda[local];
                    record.auxiliary[local] = w[local];
                }
            }
            for (int a = 0; a < 3; ++a) {
                const int row = static_cast<int>(nodes[a]);
                if (!activeNodes[static_cast<std::size_t>(row)])
                    continue;
                const Real lumpedVolume = cellArea_m2 / 3.0;
                // Eq. 231 is written for the dimensionless beta-weighted
                // driving potential w.  With C=gamma*hbar^2/(6*m*q), its
                // voltage form is div(grad(w)) + theta*|grad(w)|^2
                // + 2*Lambda/C = 0.  Lambda is already stored in volts, so
                // no additional thermal-voltage factor belongs here.
                const Real equationScale = 2.0 / cellCoefficientVm2;
                Real gradientContribution = 0.0;
                if (config.globalDiscretization == "exponential_fitted" &&
                    std::abs(config.theta - 0.5) < 1.0e-14) {
                    // For theta=1/2, div(grad(w))+|grad(w)|^2/2 =
                    // 2*laplace(u)/u with u=exp(w/2). Normalize u per cell;
                    // the arbitrary Nref then cancels without overflow.
                    const Real maximumW = std::max({w[0], w[1], w[2]});
                    Real u[3]{};
                    for (int d = 0; d < 3; ++d)
                        u[d] = safeExponential(0.5 * (w[d] - maximumW));
                    const Real inverseUa = 1.0 / u[a];
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real contribution =
                            -2.0 * stiffness * u[d] * inverseUa;
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                    }
                    Real rowWeightedU = 0.0;
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        rowWeightedU += stiffness * u[d] * inverseUa;
                    }
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        Real derivative = stiffness * u[d] * inverseUa /
                            thermalVoltage_V;
                        if (d == a)
                            derivative -= rowWeightedU / thermalVoltage_V;
                        addJacobian(
                            row, static_cast<int>(nodes[d]), derivative);
                    }
                } else if (config.globalDiscretization == "cvfem_full") {
                    // Control-volume finite element integration over the
                    // median-dual subcell.  The Laplacian is evaluated as
                    // flux through the two centroid-to-edge-midpoint faces.
                    // The nonlinear field source is constant in a Tri3 and
                    // therefore receives exactly one third of the cell area.
                    for (int d = 0; d < 3; ++d) {
                        const int col = static_cast<int>(nodes[d]);
                        const Real contribution = cvfemFlux[a][d] * w[d];
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        addJacobian(
                            row, col,
                            -cvfemFlux[a][d] / thermalVoltage_V);
                    }
                    gradientContribution =
                        config.theta * gradientSquared * lumpedVolume;
                } else {
                    for (int d = 0; d < 3; ++d) {
                        const int col = static_cast<int>(nodes[d]);
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real contribution = -stiffness * w[d];
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        addJacobian(row, col, stiffness / thermalVoltage_V);
                    }
                    gradientContribution =
                        config.theta * gradientSquared * lumpedVolume;
                }
                const Real reactionContribution =
                    equationScale * lambda[a] * lumpedVolume;
                residual(row) += gradientContribution + reactionContribution;
                if (diagnostic != nullptr) {
                    record.gradientSquared[a] = gradientContribution;
                    record.reaction[a] = reactionContribution;
                    diagnostic->gradientSquared(row) += gradientContribution;
                    diagnostic->reaction(row) += reactionContribution;
                }
                addJacobian(row, row, equationScale * lumpedVolume);
                if (config.globalDiscretization != "exponential_fitted" ||
                    std::abs(config.theta - 0.5) >= 1.0e-14) {
                    for (int d = 0; d < 3; ++d) {
                        Real gradientDerivative = 0.0;
                        for (int e = 0; e < 3; ++e) {
                            gradientDerivative += -2.0 * w[e] *
                                (b[d] * b[e] + c[d] * c[e]) /
                                (fourArea2 * thermalVoltage_V);
                        }
                        addJacobian(
                            row, static_cast<int>(nodes[d]),
                            config.theta * gradientDerivative * lumpedVolume);
                    }
                }
            }
            if (diagnostic != nullptr)
                diagnostic->cells.push_back(std::move(record));
        }
        // DEVSIM/Garcia-Asenov oxide closure.  In the oxide integrated
        // equation Lambda*V - b_n,ox*S/x_n = 0.  Multiplication by the
        // expanded Eq. 231 row scale 2/b_n,ox gives the constant face source
        // below.  Its Jacobian is exactly zero; the Lambda reaction Jacobian
        // is already assembled in the oxide cells.
        for (const auto& segment : oxideInterfaceSegments) {
            const Edge& edge = mesh.getEdge(segment.edgeId);
            const Real endpointLength_m =
                0.5 * edge.length * lengthScale_m;
            const Real boundaryContribution =
                -2.0 * endpointLength_m / oxidePenetrationDepth_m;
            const int rows[2] = {
                static_cast<int>(edge.n0), static_cast<int>(edge.n1)};
            for (const int row : rows) {
                if (!activeNodes[static_cast<std::size_t>(row)] ||
                    dirichletLambda_V.contains(static_cast<Index>(row)))
                    continue;
                residual(row) += boundaryContribution;
                if (diagnostic != nullptr)
                    diagnostic->interfaceBoundary(row) += boundaryContribution;
            }
        }
        for (int row = 0; row < nodeCount; ++row) {
            if (activeNodes[static_cast<std::size_t>(row)])
                continue;
            residual(row) = state(row);
            addJacobian(row, row, 1.0);
        }
        for (const auto& [node, lambda] : dirichletLambda_V) {
            const int row = static_cast<int>(node);
            residual(row) = state(row) -
                (lambda - nodeOutputBandDrive_V(row));
            addJacobian(row, row, 1.0);
        }
    };
    if (!config.residualDiagnosticPrefix.empty()) {
        VectorXd diagnosticResidual;
        Eq231ResidualDiagnostic diagnostic;
        assemble(continuousPotential, diagnosticResidual, nullptr, &diagnostic);
        writeEq231ResidualDiagnostic(
            config.residualDiagnosticPrefix, mesh, diagnostic,
            continuousPotential, nodeOutputBandDrive_V, activeNodes,
            dirichletLambda_V);
    }
    for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
        VectorXd residual;
        std::vector<Eigen::Triplet<Real>> triplets;
        assemble(continuousPotential, residual, &triplets, nullptr);
        SparseMatrixd jacobian(nodeCount, nodeCount);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
        for (int row = 0; row < nodeCount; ++row) {
            const auto fixed = dirichletLambda_V.find(static_cast<Index>(row));
            if (activeNodes[static_cast<std::size_t>(row)] &&
                fixed == dirichletLambda_V.end())
                continue;
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        jacobian.prune(0.0);
        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        if (initialResidualNorm < 0.0)
            initialResidualNorm = result.residualInfinityNorm;
        const Real residualTolerance = config.absoluteTolerance_V +
            config.relativeTolerance * std::max<Real>(initialResidualNorm, 1.0);
        if (result.residualInfinityNorm <= residualTolerance) {
            result.iterations = iteration;
            result.converged = true;
            break;
        }
        VectorXd update = solver.solve(jacobian, -residual);
        Real updateNorm = update.lpNorm<Eigen::Infinity>();
        Eigen::Index maxUpdateIndex = 0;
        update.cwiseAbs().maxCoeff(&maxUpdateIndex);
        if (config.maxUpdate_V > 0.0 && updateNorm > config.maxUpdate_V) {
            update *= config.maxUpdate_V / updateNorm;
            updateNorm = config.maxUpdate_V;
        }
        Real lineScale = config.damping;
        VectorXd trialResidual;
        VectorXd trial = continuousPotential + lineScale * update;
        assemble(trial, trialResidual, nullptr, nullptr);
        while (trialResidual.lpNorm<Eigen::Infinity>() >=
                   result.residualInfinityNorm && lineScale > 1.0e-6) {
            lineScale *= 0.5;
            trial = continuousPotential + lineScale * update;
            assemble(trial, trialResidual, nullptr, nullptr);
        }
        const Real trialResidualNorm =
            trialResidual.lpNorm<Eigen::Infinity>();
        // A vanishing line-search scale is a stagnation/failure condition,
        // not Newton convergence.  Previously the scaled update alone could
        // accept a state whose residual was still many orders of magnitude
        // larger than the incoming residual.
        if (!(trialResidualNorm < result.residualInfinityNorm)) {
            result.lastUpdateInfinityNorm_V = 0.0;
            result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
            result.maxUpdateNodeValue_V = 0.0;
            break;
        }
        continuousPotential = std::move(trial);
        for (const auto& [node, lambda] : dirichletLambda_V) {
            continuousPotential(static_cast<int>(node)) =
                lambda - nodeOutputBandDrive_V(static_cast<int>(node));
        }
        result.iterations = iteration;
        result.lastUpdateInfinityNorm_V = lineScale * updateNorm;
        result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
        result.maxUpdateNodeValue_V = lineScale * update(maxUpdateIndex);
        result.residualInfinityNorm = trialResidualNorm;
        const Real scale = std::max<Real>(continuousPotential.lpNorm<Eigen::Infinity>(), 1.0);
        if (lineScale * updateNorm <= config.absoluteTolerance_V +
                config.relativeTolerance * scale &&
            trialResidualNorm <= residualTolerance) {
            result.converged = true;
            break;
        }
    }
    result.potentialLike_V = continuousPotential;
    result.potential_V = continuousPotential + nodeOutputBandDrive_V;
    result.potentialInfinityNorm_V =
        result.potential_V.lpNorm<Eigen::Infinity>();
    return result;
}

} // namespace vela
