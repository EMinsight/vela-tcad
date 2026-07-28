#pragma once

/**
 * @file AssemblerUtils.h
 * @brief Shared FVM assembly helpers used by PoissonAssembler and DDAssembler.
 *
 * All functions are free (non-member) to avoid duplication between assemblers.
 * Include this header from assembler .cpp files; do not expose it as part of
 * the public library API.
 */

#include "vela/core/Types.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/ChargeSpec.h"
#include "vela/equation/Tri3LocalForwardAD.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/Material.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/BandgapNarrowing.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <string>
#include <utility>
#include <type_traits>

namespace vela::detail {

inline Real physicalPotentialCentralDifferenceStep(
    Real physicalValue, Real potentialScale, Real relativeStep)
{
    return relativeStep * std::max(potentialScale, std::abs(physicalValue));
}


// ---------------------------------------------------------------------------
// Fixed and interface charge helpers
// ---------------------------------------------------------------------------

struct RegionPairKey {
    std::string first;
    std::string second;

    bool operator==(const RegionPairKey& other) const
    {
        return first == other.first && second == other.second;
    }
};

struct RegionPairKeyHash {
    std::size_t operator()(const RegionPairKey& key) const
    {
        const std::hash<std::string> hash;
        std::size_t seed = hash(key.first);
        seed ^= hash(key.second) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

inline RegionPairKey makeRegionPairKey(std::string a, std::string b)
{
    if (b < a)
        std::swap(a, b);
    return RegionPairKey{std::move(a), std::move(b)};
}

inline Real triangleArea(const DeviceMesh& mesh, const Cell& cell)
{
    if (cell.node_ids.size() < 3) return 0.0;

    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);

    return 0.5 * std::abs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
}

inline std::unordered_map<std::string, Real> fixedChargeByRegion(
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::string& context)
{
    std::unordered_map<std::string, Real> fixedByRegion;
    for (const auto& spec : fixedCharges) {
        const auto [_, inserted] = fixedByRegion.emplace(spec.region, spec.fixedCharge);
        if (!inserted)
            throw std::invalid_argument(
                context + ": duplicate fixed_charge_m3 for region '" + spec.region + "'.");
    }
    return fixedByRegion;
}

inline std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetChargeByRegionPair(
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges)
{
    std::unordered_map<RegionPairKey, Real, RegionPairKeyHash> sheetByRegionPair;
    for (const auto& spec : sheetCharges)
        sheetByRegionPair[makeRegionPairKey(spec.region0, spec.region1)] += spec.totalSheetCharge();
    return sheetByRegionPair;
}

inline VectorXd computeFixedAndInterfaceChargeRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    const std::string& context,
    Real chargeAreaFactor = 1.0,
    Real chargeLineFactor = 1.0)
{
    VectorXd contribution = VectorXd::Zero(static_cast<int>(mesh.numNodes()));

    const auto fixedByRegion = fixedChargeByRegion(fixedCharges, context);
    if (!fixedByRegion.empty()) {
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            const Region& region = mesh.getRegion(cell.region_id);
            auto it = fixedByRegion.find(region.name);
            if (it == fixedByRegion.end()) continue;

            const Real nodeCharge = constants::q * it->second * triangleArea(mesh, cell) *
                chargeAreaFactor / 3.0;
            for (Index nid : cell.node_ids)
                contribution(static_cast<int>(nid)) += nodeCharge;
        }
    }

    const auto sheetByRegionPair = sheetChargeByRegionPair(sheetCharges);
    if (!sheetByRegionPair.empty()) {
        for (Index e = 0; e < mesh.numEdges(); ++e) {
            const auto& cells = edgeCells[e];
            if (cells.size() != 2) continue;

            const Region& r0 = mesh.getRegion(mesh.getCell(cells[0]).region_id);
            const Region& r1 = mesh.getRegion(mesh.getCell(cells[1]).region_id);
            const auto it = sheetByRegionPair.find(makeRegionPairKey(r0.name, r1.name));
            if (it == sheetByRegionPair.end()) continue;

            const Edge& edge = mesh.getEdge(e);
            const Real endpointCharge = constants::q * it->second * edge.length *
                chargeLineFactor * 0.5;
            contribution(static_cast<int>(edge.n0)) += endpointCharge;
            contribution(static_cast<int>(edge.n1)) += endpointCharge;
        }
    }

    return contribution;
}

inline void addFixedAndInterfaceChargeToRhs(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges,
    VectorXd& rhs,
    const std::string& context,
    Real chargeAreaFactor = 1.0,
    Real chargeLineFactor = 1.0)
{
    rhs += computeFixedAndInterfaceChargeRhs(
        mesh, edgeCells, fixedCharges, sheetCharges, context,
        chargeAreaFactor, chargeLineFactor);
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Return precomputed per-node control-volume areas.
inline std::vector<Real> computeNodeVolumes(const DeviceMesh& mesh)
{
    const Index N = mesh.numNodes();
    std::vector<Real> vol(N, 0.0);
    for (Index i = 0; i < N; ++i)
        vol[i] = mesh.getNode(i).volume;
    return vol;
}

/// Return precomputed per-edge box coupling lengths.
inline std::vector<Real> computeEdgeCouplings(const DeviceMesh& mesh)
{
    const Index E = mesh.numEdges();
    std::vector<Real> couple(E, 0.0);
    for (Index e = 0; e < E; ++e)
        couple[e] = mesh.getEdge(e).couple;
    return couple;
}

/// Return a per-node max adjacent-edge scalar-gradient magnitude [scalar unit/m].
inline std::vector<Real> computeNodeScalarGradientMagnitudes(const VectorXd& value,
                                                            const DeviceMesh& mesh)
{
    std::vector<Real> maxField(mesh.numNodes(), 0.0);
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        if (edge.length <= 1.0e-30)
            continue;
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real edgeField = std::abs((value(j) - value(i)) / edge.length);
        maxField[edge.n0] = std::max(maxField[edge.n0], edgeField);
        maxField[edge.n1] = std::max(maxField[edge.n1], edgeField);
    }
    return maxField;
}

/// Return a per-node Sentaurus-like cell-gradient electric-field magnitude [V/m].
inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh, Real fieldFactor = 1.0);

/// Build edge -> adjacent cell ids map.
inline std::vector<std::vector<Index>> buildEdgeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> edgeCells(mesh.numEdges());
    std::unordered_map<Index, Index> pairToEdge;
    const Index N = mesh.numNodes();
    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        pairToEdge[edge.n0 * N + edge.n1] = e;
    }
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell = mesh.getCell(c);
        if (cell.node_ids.size() < 3) continue;
        const Index nids[3] = {
            cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        for (int k = 0; k < 3; ++k) {
            Index a = nids[k];
            Index b = nids[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            auto it = pairToEdge.find(a * N + b);
            if (it != pairToEdge.end())
                edgeCells[it->second].push_back(c);
        }
    }
    return edgeCells;
}

inline std::vector<std::vector<Index>> buildCellEdgeMap(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh)
{
    std::vector<std::vector<Index>> cellEdges(static_cast<std::size_t>(mesh.numCells()));
    for (Index edgeId = 0; edgeId < edgeCells.size(); ++edgeId) {
        for (Index cellId : edgeCells[edgeId]) {
            if (cellId < mesh.numCells())
                cellEdges[static_cast<std::size_t>(cellId)].push_back(edgeId);
        }
    }
    return cellEdges;
}

inline Real cellSmoothedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               rawEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        if (cellId >= cellEdges.size())
            continue;
        Real cellSum = 0.0;
        int cellEdgeCount = 0;
        for (Index otherEdgeId : cellEdges[cellId]) {
            if (otherEdgeId >= rawEdgeFlux.size())
                continue;
            cellSum += rawEdgeFlux[otherEdgeId];
            ++cellEdgeCount;
        }
        if (cellEdgeCount <= 0)
            continue;
        edgeSum += cellSum / static_cast<Real>(cellEdgeCount);
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < rawEdgeFlux.size() ? rawEdgeFlux[edgeId] : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
}

inline Real cellVectorCurrentMagnitude(
    Index                                  cellId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (cellId >= cellEdges.size())
        return 0.0;

    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    Real absSum = 0.0;
    int used = 0;
    for (Index edgeId : cellEdges[cellId]) {
        if (edgeId >= signedEdgeFlux.size())
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.length <= 1.0e-30)
            continue;
        const Node& n0 = mesh.getNode(edge.n0);
        const Node& n1 = mesh.getNode(edge.n1);
        const Real tx = (n1.x - n0.x) / edge.length;
        const Real ty = (n1.y - n0.y) / edge.length;
        const Real flux = signedEdgeFlux[edgeId];
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
    if (used < 2 || std::abs(det) <= 1.0e-24 * scale)
        return used > 0 ? absSum / static_cast<Real>(used) : 0.0;

    const Real jx = (b0 * a11 - b1 * a01) / det;
    const Real jy = (a00 * b1 - a01 * b0) / det;
    return std::sqrt(jx * jx + jy * jy);
}

inline Real cellVectorReconstructedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        const Real cellMagnitude = cellVectorCurrentMagnitude(
            cellId, signedEdgeFlux, cellEdges, mesh);
        if (cellMagnitude <= 0.0)
            continue;
        edgeSum += cellMagnitude;
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < signedEdgeFlux.size() ? std::abs(signedEdgeFlux[edgeId]) : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
}


inline Index edgeIdForNodePair(
    const DeviceMesh&          mesh,
    const std::vector<Index>&  candidateEdges,
    Index                      a,
    Index                      b)
{
    for (Index edgeId : candidateEdges) {
        if (edgeId >= mesh.numEdges())
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        if ((edge.n0 == a && edge.n1 == b) || (edge.n0 == b && edge.n1 == a))
            return edgeId;
    }
    return mesh.numEdges();
}

inline Point2 cellCentroid(const DeviceMesh& mesh, const Cell& cell)
{
    Point2 centroid = Point2::Zero();
    if (cell.node_ids.empty())
        return centroid;
    for (Index nodeId : cell.node_ids) {
        const Node& node = mesh.getNode(nodeId);
        centroid += Point2{node.x, node.y};
    }
    return centroid / static_cast<Real>(cell.node_ids.size());
}

inline Point2 medianDualFaceNormal(
    const DeviceMesh& mesh,
    const Cell&       cell,
    Index             ownerNode,
    Index             neighborNode)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
        return Point2::Zero();
    const Node& owner = mesh.getNode(ownerNode);
    const Node& neighbor = mesh.getNode(neighborNode);
    const Point2 centroid = cellCentroid(mesh, cell);
    const Point2 midpoint{0.5 * (owner.x + neighbor.x), 0.5 * (owner.y + neighbor.y)};
    const Point2 segment = midpoint - centroid;
    const Real length = segment.norm();
    if (length <= 1.0e-30)
        return Point2::Zero();
    Point2 normal{segment.y() / length, -segment.x() / length};
    const Point2 towardNeighbor{neighbor.x - owner.x, neighbor.y - owner.y};
    if (normal.dot(towardNeighbor) < 0.0)
        normal = -normal;
    return normal;
}

inline Real medianDualFaceLength(const DeviceMesh& mesh, const Cell& cell, Index a, Index b)
{
    const Node& na = mesh.getNode(a);
    const Node& nb = mesh.getNode(b);
    const Point2 centroid = cellCentroid(mesh, cell);
    const Point2 midpoint{0.5 * (na.x + nb.x), 0.5 * (na.y + nb.y)};
    return (midpoint - centroid).norm();
}

inline Real medianDualCellVectorCurrentMagnitude(
    Index                                  cellId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (cellId >= mesh.numCells() || cellId >= cellEdges.size())
        return 0.0;
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() < 3)
        return 0.0;

    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    Real absSum = 0.0;
    int used = 0;
    for (int k = 0; k < 3; ++k) {
        const Index owner = cell.node_ids[static_cast<std::size_t>(k)];
        const Index neighbor = cell.node_ids[static_cast<std::size_t>((k + 1) % 3)];
        const Index edgeId = edgeIdForNodePair(mesh, cellEdges[cellId], owner, neighbor);
        if (edgeId >= signedEdgeFlux.size())
            continue;
        const Point2 normal = medianDualFaceNormal(mesh, cell, owner, neighbor);
        const Real normalNorm = normal.norm();
        if (normalNorm <= 1.0e-30)
            continue;
        const Edge& edge = mesh.getEdge(edgeId);
        const Real orientation = (edge.n0 == owner && edge.n1 == neighbor) ? 1.0 : -1.0;
        const Real flux = orientation * signedEdgeFlux[edgeId];
        const Real weight = std::max(medianDualFaceLength(mesh, cell, owner, neighbor), Real{1.0e-300});
        const Real nx = normal.x() / normalNorm;
        const Real ny = normal.y() / normalNorm;
        a00 += weight * nx * nx;
        a01 += weight * nx * ny;
        a11 += weight * ny * ny;
        b0 += weight * nx * flux;
        b1 += weight * ny * flux;
        absSum += std::abs(flux);
        ++used;
    }

    const Real det = a00 * a11 - a01 * a01;
    const Real scale = std::max({std::abs(a00 * a11), std::abs(a01 * a01), Real{1.0}});
    if (used < 2 || std::abs(det) <= 1.0e-24 * scale)
        return used > 0 ? absSum / static_cast<Real>(used) : 0.0;
    const Real jx = (b0 * a11 - b1 * a01) / det;
    const Real jy = (a00 * b1 - a01 * b0) / det;
    return std::sqrt(jx * jx + jy * jy);
}

inline Real medianDualFaceVectorReconstructedEdgeFluxMagnitude(
    Index                                  edgeId,
    const std::vector<Real>&               signedEdgeFlux,
    const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<std::vector<Index>>& cellEdges,
    const DeviceMesh&                      mesh)
{
    if (edgeId >= edgeCells.size())
        return 0.0;
    Real edgeSum = 0.0;
    int adjacentCellCount = 0;
    for (Index cellId : edgeCells[edgeId]) {
        const Real cellMagnitude = medianDualCellVectorCurrentMagnitude(
            cellId, signedEdgeFlux, cellEdges, mesh);
        if (cellMagnitude <= 0.0)
            continue;
        edgeSum += cellMagnitude;
        ++adjacentCellCount;
    }
    if (adjacentCellCount <= 0)
        return edgeId < signedEdgeFlux.size() ? std::abs(signedEdgeFlux[edgeId]) : 0.0;
    return edgeSum / static_cast<Real>(adjacentCellCount);
}
/// Build node -> adjacent cell ids map.
inline std::vector<std::vector<Index>> buildNodeCellMap(const DeviceMesh& mesh)
{
    std::vector<std::vector<Index>> nodeCells(mesh.numNodes());
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const Cell& cell = mesh.getCell(c);
        for (Index nodeId : cell.node_ids)
            nodeCells[nodeId].push_back(c);
    }
    return nodeCells;
}

// ---------------------------------------------------------------------------
// Material helpers
// ---------------------------------------------------------------------------

/// Return the average value of a material property over edge-adjacent cells.
/// Falls back to @p fallback only when the edge has no adjacent cells.
/// Throws if any adjacent cell references an unknown material.
inline Real edgeAvgMaterialProp(
    const std::vector<Index>& cells,
    const DeviceMesh&          mesh,
    const MaterialDatabase&    matdb,
    Real Material::*           prop,
    Real                       fallback)
{
    if (cells.empty()) return fallback;
    Real sum = 0.0;
    for (Index c : cells) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        sum += matdb.getMaterial(region.material).*prop;
    }
    return sum / static_cast<Real>(cells.size());
}


/// Build one temperature-adjusted material per mesh cell for hot-path reuse.
inline std::vector<Material> buildCellMaterials(const DeviceMesh&       mesh,
                                                const MaterialDatabase& matdb,
                                                Real                    temperature_K)
{
    std::vector<Material> materials;
    materials.reserve(mesh.numCells());
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& region = mesh.getRegion(mesh.getCell(c).region_id);
        materials.push_back(matdb.getMaterial(region.material, temperature_K));
    }
    return materials;
}


inline Real cellCentroidPotential(const DeviceMesh& mesh, const VectorXd& psi, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return std::numeric_limits<Real>::quiet_NaN();

    Real sum = 0.0;
    for (Index nodeId : cell.node_ids)
        sum += psi(static_cast<int>(nodeId));
    return sum / static_cast<Real>(cell.node_ids.size());
}

inline std::pair<Real, Real> cellCentroid(const DeviceMesh& mesh, Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return {std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};

    Real x = 0.0;
    Real y = 0.0;
    for (Index nodeId : cell.node_ids) {
        const Node& node = mesh.getNode(nodeId);
        x += node.x;
        y += node.y;
    }
    const Real invCount = 1.0 / static_cast<Real>(cell.node_ids.size());
    return {x * invCount, y * invCount};
}

inline Real estimateSurfaceNormalField(const std::vector<Index>& cells,
                                       const DeviceMesh& mesh,
                                       const VectorXd& psi,
                                       Index edgeId,
                                       Index cellId)
{
    const Edge& edge = mesh.getEdge(edgeId);
    if (edge.length <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();

    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    const Real normalX = -(n1.y - n0.y) / edge.length;
    const Real normalY =  (n1.x - n0.x) / edge.length;

    const auto [cx, cy] = cellCentroid(mesh, cellId);
    const Real cellPhi = cellCentroidPotential(mesh, psi, cellId);
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cellPhi))
        return std::numeric_limits<Real>::quiet_NaN();

    Real maxField = std::numeric_limits<Real>::quiet_NaN();
    for (Index otherCellId : cells) {
        if (otherCellId == cellId)
            continue;
        const auto [ox, oy] = cellCentroid(mesh, otherCellId);
        const Real otherPhi = cellCentroidPotential(mesh, psi, otherCellId);
        if (!std::isfinite(ox) || !std::isfinite(oy) || !std::isfinite(otherPhi))
            continue;
        const Real normalDistance = std::abs((ox - cx) * normalX + (oy - cy) * normalY);
        if (normalDistance <= 1.0e-30)
            continue;
        const Real field = std::abs((otherPhi - cellPhi) / normalDistance);
        if (!std::isfinite(maxField) || field > maxField)
            maxField = field;
    }
    if (std::isfinite(maxField))
        return maxField;

    const Real edgePhi = 0.5 * (psi(static_cast<int>(edge.n0)) + psi(static_cast<int>(edge.n1)));
    const Real mx = 0.5 * (n0.x + n1.x);
    const Real my = 0.5 * (n0.y + n1.y);
    const Real normalDistance = std::abs((cx - mx) * normalX + (cy - my) * normalY);
    if (normalDistance <= 1.0e-30)
        return std::numeric_limits<Real>::quiet_NaN();
    return std::abs((cellPhi - edgePhi) / normalDistance);
}

inline Real cellAverageTotalImpurity(const DeviceMesh& mesh,
                                     const DopingModel& doping,
                                     Index cellId)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.node_ids.empty())
        return 0.0;
    Real sum = 0.0;
    for (Index nodeId : cell.node_ids)
        sum += doping.totalImpurity(nodeId);
    return sum / static_cast<Real>(cell.node_ids.size());
}

inline Real edgeMobilityDopingConcentration(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const Edge& edge,
    Index cellId,
    const MobilityModelConfig* config)
{
    const std::string basis = config != nullptr
        ? config->dopingConcentrationBasis
        : "net_doping";
    if (basis == "total_impurity") {
        return 0.5 * (doping.totalImpurity(edge.n0) +
                      doping.totalImpurity(edge.n1));
    }
    if (basis == "cell_reconstructed_total_impurity")
        return cellAverageTotalImpurity(mesh, doping, cellId);
    return 0.5 * (doping.netDoping(edge.n0) + doping.netDoping(edge.n1));
}

inline Real nodeMobilityDopingConcentration(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    Index nodeId,
    Index cellId,
    const MobilityModelConfig* config)
{
    const std::string basis = config != nullptr
        ? config->dopingConcentrationBasis
        : "net_doping";
    if (basis == "total_impurity")
        return doping.totalImpurity(nodeId);
    if (basis == "cell_reconstructed_total_impurity")
        return cellAverageTotalImpurity(mesh, doping, cellId);
    return doping.netDoping(nodeId);
}
/// Return average model mobility [m^2/V/s] for edge @p edgeId.
inline Real edgeMobility(const std::vector<std::vector<Index>>& edgeCells,
                         const DeviceMesh&                       mesh,
                         const DopingModel&                      doping,
                         const MobilityModel&                    mobility,
                         const std::vector<Material>&            cellMaterials,
                         Index                                   edgeId,
                         CarrierType                             carrier,
                         Real                                    electricField,
                         const MobilityModelConfig*              mobilityConfig = nullptr,
                         const VectorXd*                         psi = nullptr)
{
    const auto& cells = edgeCells[edgeId];
    if (cells.empty()) return 0.0;

    const Edge& edge = mesh.getEdge(edgeId);

    const bool surfaceEnabled =
        mobilityConfig != nullptr && isSurfaceMobilityModel(*mobilityConfig);
    std::vector<std::string> adjacentRegionNames;
    if (surfaceEnabled) {
        adjacentRegionNames.reserve(cells.size());
        for (Index c : cells)
            adjacentRegionNames.push_back(mesh.getRegion(mesh.getCell(c).region_id).name);
    }

    Real sum = 0.0;
    Index contributingCells = 0;
    for (Index c : cells) {
        const Material& material = cellMaterials.at(static_cast<std::size_t>(c));
        const Real baseMobility = (carrier == CarrierType::Electron) ? material.mun : material.mup;
        if (baseMobility <= 0.0)
            continue;

        // Average only transport-capable cells.  This keeps oxide-only
        // edges pinned while preserving lateral semiconductor transport on
        // edges that lie along a semiconductor/oxide interface. Surface
        // mobility is enabled only on configured regions/interfaces; when no
        // normal-field estimate is available the NaN field disables the surface
        // factor while preserving any high-field velocity saturation.
        const Region& region = mesh.getRegion(mesh.getCell(c).region_id);
        const bool surfaceApplies = surfaceEnabled &&
            surfaceMobilityAppliesToRegionPair(*mobilityConfig, region.name, adjacentRegionNames);
        const Real surfaceNormalField = (surfaceApplies && psi != nullptr)
            ? estimateSurfaceNormalField(cells, mesh, *psi, edgeId, c)
            : std::numeric_limits<Real>::quiet_NaN();
        const Real mobilityDoping = edgeMobilityDopingConcentration(
            mesh, doping, edge, c, mobilityConfig);
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, mobilityDoping, 0.0, 0.0, electricField, surfaceNormalField)
            : mobility.holeMobility(
                material, mobilityDoping, 0.0, 0.0, electricField, surfaceNormalField);
        if (modelMobility <= 0.0)
            continue;

        sum += modelMobility;
        ++contributingCells;
    }
    if (contributingCells == 0)
        return 0.0;
    return sum / static_cast<Real>(contributingCells);
}

/// Return average model mobility [m^2/V/s] over semiconductor cells adjacent to a node.
inline Real nodeMobility(const std::vector<std::vector<Index>>& nodeCells,
                         const DeviceMesh&                     mesh,
                         const DopingModel&                    doping,
                         const MobilityModel&                  mobility,
                         const std::vector<Material>&          cellMaterials,
                         Index                                 nodeId,
                         CarrierType                           carrier,
                         Real                                  drivingField,
                         const MobilityModelConfig*            mobilityConfig = nullptr)
{
    const auto& cells = nodeCells[nodeId];
    if (cells.empty())
        return 0.0;

    Real sum = 0.0;
    Index contributingCells = 0;
    for (Index c : cells) {
        const Material& material = cellMaterials.at(static_cast<std::size_t>(c));
        const Real baseMobility = (carrier == CarrierType::Electron) ? material.mun : material.mup;
        if (baseMobility <= 0.0)
            continue;
        const Real mobilityDoping = nodeMobilityDopingConcentration(
            mesh, doping, nodeId, c, mobilityConfig);
        const Real modelMobility = (carrier == CarrierType::Electron)
            ? mobility.electronMobility(
                material, mobilityDoping, 0.0, 0.0, drivingField)
            : mobility.holeMobility(
                material, mobilityDoping, 0.0, 0.0, drivingField);
        if (modelMobility <= 0.0)
            continue;
        sum += modelMobility;
        ++contributingCells;
    }
    if (contributingCells == 0)
        return 0.0;
    return sum / static_cast<Real>(contributingCells);
}

inline Real interpolatedAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                              Real                               drivingField,
                                              Real                               electricField,
                                              Real                               carrierDensity,
                                              Real                               referenceDensity)
{
    if (config.debugRawVanOverstraeten ||
        config.drivingForceInterpolation != "quasi_fermi_to_electric_field" ||
        referenceDensity <= 0.0) {
        return drivingField;
    }
    const Real carrier = std::max(carrierDensity, 0.0);
    const Real weight = carrier / (carrier + referenceDensity);
    return weight * drivingField + (1.0 - weight) * electricField;
}

inline Real electronAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                          Real                               drivingField,
                                          Real                               electricField,
                                          Real                               electronDensity)
{
    return interpolatedAvalancheDrivingField(
        config,
        drivingField,
        electricField,
        electronDensity,
        config.electronDrivingForceRefDensity);
}

inline Real holeAvalancheDrivingField(const ImpactIonizationModelConfig& config,
                                      Real                               drivingField,
                                      Real                               electricField,
                                      Real                               holeDensity)
{
    return interpolatedAvalancheDrivingField(
        config,
        drivingField,
        electricField,
        holeDensity,
        config.holeDrivingForceRefDensity);
}

inline bool usesCurrentAlignedAvalancheDrivingForce(
    const ImpactIonizationModelConfig& config)
{
    return config.drivingForce == "grad_potential_parallel_j" ||
           config.drivingForce == "effective_field_parallel_j";
}

inline bool usesQuasiFermiAvalancheDrivingForce(
    const ImpactIonizationModelConfig& config)
{
    return config.debugRawVanOverstraeten ||
           config.drivingForce == "quasi_fermi_gradient";
}

inline Real parallelCurrentAvalancheDrivingField(Real signedDrivingField,
                                                Real signedCurrentProxy)
{
    if (!std::isfinite(signedDrivingField) || !std::isfinite(signedCurrentProxy) ||
        std::abs(signedCurrentProxy) <= 0.0) {
        return 0.0;
    }
    const Real currentSign = signedCurrentProxy > 0.0 ? 1.0 : -1.0;
    return std::max(signedDrivingField * currentSign, 0.0);
}

/// Resolves the legacy SG edge-current avalanche source-volume factor used in
/// `factor * h * edge.couple`. A finite `source_volume_factor` overrides the
/// named `source_volume_policy` preset; `0` falls back to the preset.
inline Real avalancheSourceVolumeFactor(const ImpactIonizationModelConfig& config)
{
    if (config.sourceVolumeFactor > 0.0)
        return config.sourceVolumeFactor;
    return config.sourceVolumePolicy == "edge_box" ? 1.0 : 0.5;
}

inline Real triangleSignedDoubleArea(const Point2& a, const Point2& b, const Point2& c)
{
    return (b.x() - a.x()) * (c.y() - a.y()) -
           (c.x() - a.x()) * (b.y() - a.y());
}

inline Point2 meshPoint(const DeviceMesh& mesh, Index node)
{
    const Node& n = mesh.getNode(node);
    return Point2{n.x, n.y};
}

inline int tri3LocalEdgeIndex(const Cell& cell, Index edgeNode0, Index edgeNode1)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return -1;
    for (int local = 0; local < 3; ++local) {
        const Index a = cell.node_ids[static_cast<std::size_t>(local)];
        const Index b = cell.node_ids[static_cast<std::size_t>((local + 1) % 3)];
        if ((a == edgeNode0 && b == edgeNode1) ||
            (a == edgeNode1 && b == edgeNode0)) {
            return local;
        }
    }
    return -1;
}

inline Real angleBetween(const Point2& a, const Point2& b)
{
    const Real denom = a.norm() * b.norm();
    if (denom <= 1.0e-300)
        return 0.0;
    const Real cosTheta = std::clamp(a.dot(b) / denom, -1.0, 1.0);
    return std::acos(cosTheta);
}

inline Real geniusTri3TruncatedPartialVolumeWithEdge(
    const DeviceMesh& mesh,
    const Cell&       cell,
    Index             edgeNode0,
    Index             edgeNode1)
{
    const int localEdge = tri3LocalEdgeIndex(cell, edgeNode0, edgeNode1);
    if (localEdge < 0)
        return 0.0;

    const std::array<Index, 3> ids = {
        cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
    const std::array<Point2, 3> p = {
        meshPoint(mesh, ids[0]), meshPoint(mesh, ids[1]), meshPoint(mesh, ids[2])};
    const Real lengthScale = std::max({
        (p[0] - p[1]).norm(), (p[1] - p[2]).norm(), (p[2] - p[0]).norm()});
    const Real distanceTolerance =
        64.0 * std::numeric_limits<Real>::epsilon() * lengthScale;
    const Real det = triangleSignedDoubleArea(p[0], p[1], p[2]);
    if (std::abs(det) <= 1.0e-300)
        return 0.0;

    const Real a2 = p[0].squaredNorm();
    const Real b2 = p[1].squaredNorm();
    const Real c2 = p[2].squaredNorm();
    const Real invDenom = 1.0 / (2.0 * det);
    const Point2 circumcenter{
        (a2 * (p[1].y() - p[2].y()) +
         b2 * (p[2].y() - p[0].y()) +
         c2 * (p[0].y() - p[1].y())) * invDenom,
        (a2 * (p[2].x() - p[1].x()) +
         b2 * (p[0].x() - p[2].x()) +
         c2 * (p[1].x() - p[0].x())) * invDenom};

    constexpr int sideNodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    std::array<Real, 3> lengths = {0.0, 0.0, 0.0};
    std::array<Real, 3> dt = {0.0, 0.0, 0.0};
    int obtuseEdge = -1;
    for (int local = 0; local < 3; ++local) {
        const Point2& p1 = p[sideNodes[local][0]];
        const Point2& p2 = p[sideNodes[local][1]];
        const Point2& p3 = p[(2 + local) % 3];
        const Point2 sideCenter = 0.5 * (p1 + p2);
        lengths[static_cast<std::size_t>(local)] = (p1 - p2).norm();
        const Real rawDistance = (sideCenter - circumcenter).norm();
        const Real distance =
            rawDistance <= distanceTolerance ? 0.0 : rawDistance;
        if ((p1 - p3).dot(p2 - p3) < 0.0) {
            dt[static_cast<std::size_t>(local)] = -distance;
            obtuseEdge = local;
        } else {
            dt[static_cast<std::size_t>(local)] = distance;
        }
    }

    if (obtuseEdge >= 0) {
        const int obtuseNode = (2 + obtuseEdge) % 3;
        const Point2& p1 = p[sideNodes[obtuseEdge][0]];
        const Point2& p2 = p[sideNodes[obtuseEdge][1]];
        const Point2& p3 = p[obtuseNode];
        const Real theta1 = angleBetween(p2 - p1, p3 - p1);
        const Real theta2 = angleBetween(p1 - p2, p3 - p2);
        const Real cos1 = std::cos(theta1);
        const Real cos2 = std::cos(theta2);
        const Point2 preEdgeCenter = 0.5 * (p1 + p3);
        const Point2 posEdgeCenter = 0.5 * (p2 + p3);
        dt[static_cast<std::size_t>(obtuseEdge)] = 0.0;
        if (std::abs(cos1) > 1.0e-300 && (p2 - p1).norm() > 1.0e-300) {
            const Point2 m1 = p1 + (p2 - p1).normalized() *
                ((preEdgeCenter - p1).norm() / cos1);
            dt[static_cast<std::size_t>((obtuseEdge + 2) % 3)] =
                (preEdgeCenter - m1).norm();
        }
        if (std::abs(cos2) > 1.0e-300 && (p1 - p2).norm() > 1.0e-300) {
            const Point2 m2 = p2 + (p1 - p2).normalized() *
                ((posEdgeCenter - p2).norm() / cos2);
            dt[static_cast<std::size_t>((obtuseEdge + 1) % 3)] =
                (posEdgeCenter - m2).norm();
        }
    }

    return 0.5 * lengths[static_cast<std::size_t>(localEdge)] *
           std::max(0.0, dt[static_cast<std::size_t>(localEdge)]);
}

inline std::array<Real, 3> tri3ElementEdgeBoxPartialVolumes(
    const DeviceMesh& mesh,
    const Cell&       cell)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux reconstruction requires a Tri3 cell");
    }
    std::array<Real, 3> partialVolumes{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 =
            cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        partialVolumes[static_cast<std::size_t>(localEdge)] =
            geniusTri3TruncatedPartialVolumeWithEdge(
                mesh, cell, node0, node1);
    }
    const Real exactArea = 0.5 * std::abs(triangleSignedDoubleArea(
        meshPoint(mesh, cell.node_ids[0]),
        meshPoint(mesh, cell.node_ids[1]),
        meshPoint(mesh, cell.node_ids[2])));
    const Real truncatedArea =
        partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
    const Real closureTolerance =
        128.0 * std::numeric_limits<Real>::epsilon() * exactArea;
    if (exactArea > 0.0 && truncatedArea > 0.0 &&
        std::abs(truncatedArea - exactArea) > closureTolerance) {
        const Real conservativeScale = exactArea / truncatedArea;
        for (Real& partialVolume : partialVolumes)
            partialVolume *= conservativeScale;
    }
    return partialVolumes;
}

inline std::array<Real, 3> tri3ElementVertexBoxMeasures(
    const DeviceMesh& mesh,
    const Cell&       cell)
{
    const auto edgePartialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    std::array<Real, 3> vertexMeasures{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Real halfMeasure =
            0.5 * edgePartialVolumes[static_cast<std::size_t>(localEdge)];
        vertexMeasures[static_cast<std::size_t>(localEdge)] += halfMeasure;
        vertexMeasures[static_cast<std::size_t>((localEdge + 1) % 3)] +=
            halfMeasure;
    }
    return vertexMeasures;
}

inline Point2 solveTri3EdgeProjectionPair(
    const Point2& tangentA,
    Real          valueA,
    const Point2& tangentB,
    Real          valueB)
{
    const Real determinant =
        tangentA.x() * tangentB.y() - tangentA.y() * tangentB.x();
    if (std::abs(determinant) <= 1.0e-14) {
        throw std::invalid_argument(
            "parallel triangle edges cannot reconstruct a current vector");
    }
    return Point2{
        (valueA * tangentB.y() - tangentA.y() * valueB) / determinant,
        (tangentA.x() * valueB - valueA * tangentB.x()) / determinant};
}

inline Point2 gssLauxTri3CurrentVector(
    const DeviceMesh&          mesh,
    const Cell&                cell,
    const std::array<Real, 3>& signedEdgeCurrent)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux reconstruction requires a Tri3 cell");
    }

    std::array<Point2, 3> tangents{};
    std::array<Real, 3> lengths{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Point2 delta =
            meshPoint(
                mesh,
                cell.node_ids[static_cast<std::size_t>(
                    (localEdge + 1) % 3)]) -
            meshPoint(
                mesh,
                cell.node_ids[static_cast<std::size_t>(localEdge)]);
        lengths[static_cast<std::size_t>(localEdge)] = delta.norm();
        if (lengths[static_cast<std::size_t>(localEdge)] <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot reconstruct a current vector");
        }
        tangents[static_cast<std::size_t>(localEdge)] =
            delta / lengths[static_cast<std::size_t>(localEdge)];
    }

    const std::array<Point2, 3> pairVectors = {
        solveTri3EdgeProjectionPair(
            tangents[0], signedEdgeCurrent[0],
            tangents[1], signedEdgeCurrent[1]),
        solveTri3EdgeProjectionPair(
            tangents[0], signedEdgeCurrent[0],
            tangents[2], signedEdgeCurrent[2]),
        solveTri3EdgeProjectionPair(
            tangents[1], signedEdgeCurrent[1],
            tangents[2], signedEdgeCurrent[2])};
    const auto pairIndex = [](int first, int second) {
        const int low = std::min(first, second);
        const int high = std::max(first, second);
        return low == 0 ? high - 1 : 2;
    };

    const auto partialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    std::array<Point2, 3> edgeVectors{};
    for (int target = 0; target < 3; ++target) {
        std::array<int, 2> others{};
        int next = 0;
        for (int candidate = 0; candidate < 3; ++candidate) {
            if (candidate != target)
                others[static_cast<std::size_t>(next++)] = candidate;
        }
        const int first = others[0];
        const int second = others[1];
        const Real firstWeight =
            2.0 * partialVolumes[static_cast<std::size_t>(first)] /
            lengths[static_cast<std::size_t>(first)];
        const Real secondWeight =
            2.0 * partialVolumes[static_cast<std::size_t>(second)] /
            lengths[static_cast<std::size_t>(second)];
        const Point2& firstPair =
            pairVectors[static_cast<std::size_t>(pairIndex(target, first))];
        const Point2& secondPair =
            pairVectors[static_cast<std::size_t>(pairIndex(target, second))];
        const Real weightSum = firstWeight + secondWeight;
        if (weightSum > 0.0) {
            edgeVectors[static_cast<std::size_t>(target)] =
                (firstWeight * firstPair + secondWeight * secondPair) /
                weightSum;
        } else {
            edgeVectors[static_cast<std::size_t>(target)] =
                0.5 * (firstPair + secondPair);
        }
    }

    const Real totalVolume =
        partialVolumes[0] + partialVolumes[1] + partialVolumes[2];
    if (totalVolume <= 1.0e-300) {
        throw std::invalid_argument(
            "triangle has no positive box partial volume");
    }
    return (
        partialVolumes[0] * edgeVectors[0] +
        partialVolumes[1] * edgeVectors[1] +
        partialVolumes[2] * edgeVectors[2]) / totalVolume;
}

inline Real geniusTruncatedEdgeSourceVolume(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId)
{
    if (edgeId >= mesh.numEdges() || edgeId >= edgeCells.size())
        return 0.0;
    const Edge& edge = mesh.getEdge(edgeId);
    Real volume = 0.0;
    for (Index cellId : edgeCells[edgeId]) {
        if (cellId >= mesh.numCells())
            continue;
        volume += geniusTri3TruncatedPartialVolumeWithEdge(
            mesh, mesh.getCell(cellId), edge.n0, edge.n1);
    }
    return volume;
}

inline Real avalancheSourceEdgeArea(
    const ImpactIonizationModelConfig&     config,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId)
{
    if (edgeId >= mesh.numEdges())
        return 0.0;
    const Edge& edge = mesh.getEdge(edgeId);
    Real area = 0.0;
    if (config.sourceVolumeFactor > 0.0) {
        area = config.sourceVolumeFactor * edge.length * edge.couple;
    } else if (config.sourceVolumePolicy == "genius_truncated") {
        area = geniusTruncatedEdgeSourceVolume(edgeCells, mesh, edgeId);
    } else {
        area = avalancheSourceVolumeFactor(config) * edge.length * edge.couple;
    }
    return area * config.sourceGeometryScale;
}

/// Validates the impact-ionization configuration shared by the Gummel and
/// Newton solver config loaders. `context` is prefixed to any thrown message.
inline void validateImpactIonizationDrivingForce(const ImpactIonizationModelConfig& config,
                                                 const char* context)
{
    const bool configuredCurrentAlignedDrivingForce =
        usesCurrentAlignedAvalancheDrivingForce(config);
    const bool currentAlignedDrivingForce =
        !config.debugRawVanOverstraeten && configuredCurrentAlignedDrivingForce;
    if (config.contactElectricFieldFallbackScope != "contact_node_cell" &&
        config.contactElectricFieldFallbackScope != "contact_boundary_face") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback_scope must be "
            "'contact_node_cell' or 'contact_boundary_face'.");
    }

    const bool validContactFallbackMode =
        config.contactElectricFieldFallbackMode == "cell_gradient_magnitude" ||
        config.contactElectricFieldFallbackMode == "face_normal" ||
        config.contactElectricFieldFallbackMode == "one_sided" ||
        config.contactElectricFieldFallbackMode == "distance_weighted_blend";
    if (!validContactFallbackMode) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback_mode must be "
            "'cell_gradient_magnitude', 'face_normal', 'one_sided', or "
            "'distance_weighted_blend'.");
    }
    if (config.contactElectricFieldFallbackMode != "cell_gradient_magnitude" &&
        config.contactElectricFieldFallbackScope != "contact_boundary_face") {
        throw std::invalid_argument(
            std::string(context) +
            ": contact-face fallback modes require scope 'contact_boundary_face'.");
    }

    if (config.contactElectricFieldFallback &&
        (config.drivingForce != "quasi_fermi_gradient" ||
         config.quasiFermiGradientDiscretization != "cell_gradient" ||
         config.generation != "current_density" ||
         (config.sourceMappingMode != "triangle_gss_gradqf_truncated" &&
          config.sourceMappingMode != "element_vertex_box_measure"))) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.contact_electric_field_fallback requires "
            "current-density GradQF with cell_gradient and a triangle or "
            "element-vertex-box source mapping.");
    }

    if (config.drivingForce != "electric_field" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !configuredCurrentAlignedDrivingForce) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force must be 'electric_field', "
            "'quasi_fermi_gradient', 'grad_potential_parallel_j', or "
            "'effective_field_parallel_j'.");
    }
    if (config.generation != "carrier_density" &&
        config.generation != "current_density") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.generation must be 'carrier_density' or "
            "'current_density'.");
    }
    if (config.currentApproximation != "mobility_density_gradient" &&
        config.currentApproximation != "density_gradient" &&
        config.currentApproximation != "grad_qf" &&
        config.currentApproximation != "cell_reconstructed" &&
        config.currentApproximation != "psi_gradient_proxy" &&
        config.currentApproximation != "cell_current_reconstructed" &&
        config.currentApproximation != "cell_vector_current_reconstructed" &&
        config.currentApproximation != "element_edge_sg_gss_laux" &&
        config.currentApproximation != "conserved_total_current") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.current_approximation must be "
            "'mobility_density_gradient', 'density_gradient', 'grad_qf', "
            "'cell_reconstructed', 'psi_gradient_proxy', 'cell_current_reconstructed', "
            "'cell_vector_current_reconstructed', 'element_edge_sg_gss_laux', "
            "or 'conserved_total_current'.");
    }
    if (config.currentMagnitudeMode != "edge_scalar_abs" &&
        config.currentMagnitudeMode != "dual_face_vector_mag") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.current_magnitude_mode must be "
            "'edge_scalar_abs' or 'dual_face_vector_mag'.");
    }
    if (config.cellReconstructedMidpointDensity != "bernoulli" &&
        config.cellReconstructedMidpointDensity != "arithmetic" &&
        config.cellReconstructedMidpointDensity != "gss_logistic") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.cell_reconstructed_midpoint_density must be "
            "'bernoulli', 'arithmetic', or 'gss_logistic'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForceInterpolation != "quasi_fermi_to_electric_field") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation.mode must be "
            "'none' or 'quasi_fermi_to_electric_field'.");
    }
    if (currentAlignedDrivingForce &&
        (config.generation != "current_density" ||
         (config.currentApproximation != "density_gradient" &&
          config.currentApproximation != "grad_qf"))) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization current-aligned driving forces require "
            "generation='current_density' with current_approximation='density_gradient' "
            "or 'grad_qf'.");
    }
    if (config.drivingForceInterpolation != "none" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !config.debugRawVanOverstraeten) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation requires "
            "driving_force='quasi_fermi_gradient'.");
    }
    if (config.quasiFermiGradientDiscretization != "edge_difference" &&
        config.quasiFermiGradientDiscretization != "cell_gradient") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_gradient_discretization must be "
            "'edge_difference' or 'cell_gradient'.");
    }
    if (config.quasiFermiGradientDiscretization == "cell_gradient" &&
        config.drivingForce != "quasi_fermi_gradient" &&
        !config.debugRawVanOverstraeten) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_gradient_discretization='cell_gradient' "
            "requires driving_force='quasi_fermi_gradient'.");
    }
    if (!std::isfinite(config.electronDrivingForceRefDensity) ||
        !std::isfinite(config.holeDrivingForceRefDensity) ||
        config.electronDrivingForceRefDensity < 0.0 ||
        config.holeDrivingForceRefDensity < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization driving-force reference densities must be "
            "finite and non-negative.");
    }
    if (!std::isfinite(config.sourceGeometryScale) ||
        config.sourceGeometryScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_geometry_scale must be positive and finite.");
    }
    if (config.sourceVolumePolicy != "genius_truncated" &&
        config.sourceVolumePolicy != "edge_half_box" &&
        config.sourceVolumePolicy != "edge_box") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_policy must be 'genius_truncated', 'edge_half_box', or 'edge_box'.");
    }
    if (config.sourceVolumeFactor != 0.0 &&
        (!std::isfinite(config.sourceVolumeFactor) ||
         config.sourceVolumeFactor < 0.5 ||
         config.sourceVolumeFactor > 1.0)) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_volume_factor must be 0 or within [0.5, 1.0].");
    }
    if (config.sourceMappingMode != "node_F_node_alpha_node_G" &&
        config.sourceMappingMode != "edge_F_edge_alpha_edge_G_to_node" &&
        config.sourceMappingMode != "cell_F_cell_alpha_cell_G_to_node" &&
        config.sourceMappingMode != "triangle_gss_gradqf_truncated" &&
        config.sourceMappingMode != "element_vertex_box_measure") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.source_mapping_mode must be 'node_F_node_alpha_node_G', "
            "'edge_F_edge_alpha_edge_G_to_node', 'cell_F_cell_alpha_cell_G_to_node', "
            "'triangle_gss_gradqf_truncated', or 'element_vertex_box_measure'.");
    }
    const bool elementEdgeGssLauxRequested =
        config.currentApproximation == "element_edge_sg_gss_laux" ||
        config.sourceMappingMode == "element_vertex_box_measure";
    if (elementEdgeGssLauxRequested &&
        (config.generation != "current_density" ||
         (config.drivingForce != "quasi_fermi_gradient" &&
          config.drivingForce != "electric_field") ||
         config.currentApproximation != "element_edge_sg_gss_laux" ||
         (config.drivingForce == "quasi_fermi_gradient" &&
          config.quasiFermiGradientDiscretization != "cell_gradient") ||
         config.sourceMappingMode != "element_vertex_box_measure")) {
        throw std::invalid_argument(
            std::string(context) +
            ": element_edge_sg_gss_laux requires the canonical element-box configuration.");
    }
    if (config.sourceMappingMode == "triangle_gss_gradqf_truncated" &&
        (config.generation != "current_density" ||
         config.drivingForce != "quasi_fermi_gradient" ||
         config.currentApproximation != "cell_reconstructed" ||
         config.currentMagnitudeMode != "edge_scalar_abs" ||
         config.cellReconstructedMidpointDensity != "gss_logistic" ||
         config.quasiFermiGradientDiscretization != "cell_gradient" ||
         config.sourceVolumePolicy != "genius_truncated" ||
         config.sourceVolumeFactor != 0.0 ||
         config.sourceGeometryScale != 1.0 ||
         config.edgeSourcePartition != "symmetric" ||
         config.drivingForceInterpolation != "none" ||
         config.quasiFermiCarrierTruncation != 0.0 ||
         config.minimumField != 0.0 ||
         config.electronDrivingForceRefDensity != 0.0 ||
         config.holeDrivingForceRefDensity != 0.0)) {
        throw std::invalid_argument(
            std::string(context) +
            ": triangle_gss_gradqf_truncated requires the canonical GSS GradQf configuration.");
    }
    if (config.cellReconstructedMidpointDensity == "gss_logistic" &&
        config.sourceMappingMode != "triangle_gss_gradqf_truncated") {
        throw std::invalid_argument(
            std::string(context) +
            ": cell_reconstructed_midpoint_density='gss_logistic' requires "
            "source_mapping_mode='triangle_gss_gradqf_truncated'.");
    }
    if (config.edgeSourcePartition != "symmetric" &&
        config.edgeSourcePartition != "qf_gradient") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.edge_source_partition must be 'symmetric' or 'qf_gradient'.");
    }
    if (!std::isfinite(config.quasiFermiCarrierTruncation) ||
        config.quasiFermiCarrierTruncation < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.quasi_fermi_carrier_truncation must be non-negative and finite.");
    }
    if (!std::isfinite(config.minimumField) || config.minimumField < 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.minimum_field_V_m must be non-negative and finite.");
    }
    if (!std::isfinite(config.aScale) || config.aScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.A_scale must be positive and finite.");
    }
    if (config.aScale != 1.0 && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.A_scale requires model='van_overstraeten'.");
    }
    if (!std::isfinite(config.bScale) || config.bScale <= 0.0) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.B_scale must be positive and finite.");
    }
    if (config.bScale != 1.0 && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.B_scale requires model='van_overstraeten'.");
    }
    if (config.debugRawVanOverstraeten && config.model != "van_overstraeten") {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.debug_raw_vanoverstraeten requires model='van_overstraeten'.");
    }
}

inline bool usesTriangleGssAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.sourceMappingMode == "triangle_gss_gradqf_truncated";
}

inline bool usesElementEdgeGssLauxAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "element_edge_sg_gss_laux" &&
           config.sourceMappingMode == "element_vertex_box_measure";
}

inline bool usesDensityGradientAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "density_gradient";
}

inline bool usesCellReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_reconstructed";
}

inline bool usesPsiGradientProxyAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "psi_gradient_proxy";
}

inline bool usesCellCurrentReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_current_reconstructed";
}

inline bool usesCellVectorCurrentReconstructedAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "cell_vector_current_reconstructed";
}

/// Avalanche source driven by the conserved total-current magnitude |F_p-F_n|
/// on each edge instead of the per-carrier local-density SG flux. The total
/// charge current is divergence-free in the converged state, so it does not
/// collapse on the depleted side of a reverse-biased junction where the
/// per-carrier flux (and hence the generation seed) otherwise vanishes.
inline bool usesConservedTotalCurrentAvalancheCurrent(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           config.currentApproximation == "conserved_total_current";
}

inline Real conservedTotalCurrentFluxMagnitude(Real electronContinuityFlux,
                                               Real holeContinuityFlux)
{
    return std::abs(holeContinuityFlux - electronContinuityFlux);
}

inline Real reconstructedAvalancheCurrentDensityMagnitude(Real mobility,
                                                         Real carrierDensity,
                                                         Real drivingField)
{
    if (mobility <= 0.0)
        return 0.0;
    return mobility * std::max(carrierDensity, 0.0) * std::abs(drivingField);
}

inline Real selectAvalancheCurrentFluxProxy(
    const ImpactIonizationModelConfig& config,
    Real rawFluxMagnitude,
    Real reconstructedFluxMagnitude,
    Real mobility,
    Real midpointDensity,
    Real impactField,
    Real electricField,
    Real conservedTotalFluxMagnitude)
{
    const bool reconstructedCurrent =
        config.currentMagnitudeMode == "dual_face_vector_mag" ||
        usesCellCurrentReconstructedAvalancheCurrent(config) ||
        usesCellVectorCurrentReconstructedAvalancheCurrent(config);
    if (reconstructedCurrent)
        return reconstructedFluxMagnitude;
    if (usesPsiGradientProxyAvalancheCurrent(config)) {
        return reconstructedAvalancheCurrentDensityMagnitude(
            mobility, midpointDensity, electricField);
    }
    if (usesCellReconstructedAvalancheCurrent(config)) {
        return reconstructedAvalancheCurrentDensityMagnitude(
            mobility, midpointDensity, impactField);
    }
    if (usesConservedTotalCurrentAvalancheCurrent(config))
        return conservedTotalFluxMagnitude;
    return rawFluxMagnitude;
}

/// Fermi/logistic weight aux2(x) = 1 / (1 + exp(x)); numerically stable.
inline Real avalancheMidpointAux2(Real x)
{
    if (x >= 0.0)
        return 1.0 / (1.0 + std::exp(x));
    const Real ex = std::exp(x);
    return 1.0 / (1.0 + ex);
}

/// Bernoulli/exponentially weighted edge-midpoint carrier density:
///   n_mid = n_i * aux2((V_i - V_j) / (2 Vt)) + n_j * aux2((V_j - V_i) / (2 Vt))
/// with aux2(x) = 1/(1+exp(x)); the two weights sum to 1. For electrons pass the
/// electrostatic potentials as (V_i, V_j); for holes swap them so the potential
/// enters with the opposite sign.
inline Real bernoulliWeightedMidpointDensity(Real density_i,
                                             Real density_j,
                                             Real potential_i,
                                             Real potential_j,
                                             Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_i - potential_j) / (2.0 * Vt);
    const Real weight_i = avalancheMidpointAux2(arg);
    const Real weight_j = avalancheMidpointAux2(-arg);
    return density_i * weight_i + density_j * weight_j;
}

inline Real cellReconstructedAvalancheMidpointDensity(
    const ImpactIonizationModelConfig& config,
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (config.cellReconstructedMidpointDensity == "arithmetic")
        return 0.5 * (density_i + density_j);
    return bernoulliWeightedMidpointDensity(density_i, density_j, potential_i, potential_j, Vt);
}

inline Real gssElectronAvalancheMidpointDensity(
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_j - potential_i) / (2.0 * Vt);
    return density_i * avalancheMidpointAux2(arg) +
           density_j * avalancheMidpointAux2(-arg);
}

inline Real gssHoleAvalancheMidpointDensity(
    Real density_i,
    Real density_j,
    Real potential_i,
    Real potential_j,
    Real Vt)
{
    if (Vt <= 0.0)
        return 0.5 * (density_i + density_j);
    const Real arg = (potential_i - potential_j) / (2.0 * Vt);
    return density_i * avalancheMidpointAux2(arg) +
           density_j * avalancheMidpointAux2(-arg);
}

inline bool usesEdgeCurrentAvalancheSource(
    const ImpactIonizationModelConfig& config)
{
    return config.generation == "current_density" &&
           (config.currentApproximation == "density_gradient" ||
            config.currentApproximation == "grad_qf" ||
            config.currentApproximation == "cell_reconstructed" ||
            config.currentApproximation == "psi_gradient_proxy" ||
            config.currentApproximation == "cell_current_reconstructed" ||
            config.currentApproximation == "cell_vector_current_reconstructed" ||
            config.currentApproximation == "element_edge_sg_gss_laux" ||
            config.currentApproximation == "conserved_total_current");
}

inline bool usesDirectionalEdgeAvalancheSourcePartition(
    const ImpactIonizationModelConfig& config)
{
    return usesEdgeCurrentAvalancheSource(config) &&
           (config.currentApproximation == "grad_qf" ||
            config.edgeSourcePartition == "qf_gradient");
}

inline bool usesQuasiFermiCarrierTruncation(const ImpactIonizationModelConfig& config)
{
    return !config.debugRawVanOverstraeten &&
           config.quasiFermiCarrierTruncation > 0.0;
}

inline bool usesCellGradientQuasiFermiAvalancheDrive(
    const ImpactIonizationModelConfig& config)
{
    return usesQuasiFermiAvalancheDrivingForce(config) &&
           config.quasiFermiGradientDiscretization == "cell_gradient";
}

inline std::vector<bool> contactNodeMask(const DeviceMesh& mesh)
{
    std::vector<bool> mask(static_cast<std::size_t>(mesh.numNodes()), false);
    for (const Contact& contact : mesh.contacts()) {
        for (Index nodeId : contact.node_ids) {
            if (nodeId < mesh.numNodes())
                mask[static_cast<std::size_t>(nodeId)] = true;
        }
    }
    return mask;
}
inline bool cellTouchesContact(const DeviceMesh& mesh, const Cell& cell)
{
    for (const Contact& contact : mesh.contacts()) {
        for (Index contactNode : contact.node_ids) {
            for (Index cellNode : cell.node_ids) {
                if (contactNode == cellNode)
                    return true;
            }
        }
    }
    return false;
}

inline bool cellTouchesContactBoundaryFace(
    const DeviceMesh& mesh,
    const Cell& cell)
{
    for (const Contact& contact : mesh.contacts()) {
        std::size_t contactVertexCount = 0;
        for (Index cellNode : cell.node_ids) {
            for (Index contactNode : contact.node_ids) {
                if (cellNode == contactNode) {
                    ++contactVertexCount;
                    break;
                }
            }
        }
        if (contactVertexCount >= 2)
            return true;
    }
    return false;
}

struct ContactBoundaryFaceGeometry {
    std::array<std::size_t, 2> faceLocalNodes{};
    std::size_t oppositeLocalNode = 0;
    Point2 unitNormal = Point2::Zero();
    Point2 faceMidpoint = Point2::Zero();
    Real oppositeToFaceMidpointDistance = 0.0;
    Real normalHeight = 0.0;
    Real centroidElectricWeight = 0.0;
};

inline bool contactBoundaryFaceGeometry(
    const DeviceMesh& mesh,
    const Cell& cell,
    ContactBoundaryFaceGeometry& geometry)
{
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return false;
    for (const Contact& contact : mesh.contacts()) {
        std::array<bool, 3> onContact{};
        for (std::size_t local = 0; local < 3; ++local) {
            onContact[local] = std::find(
                contact.node_ids.begin(), contact.node_ids.end(),
                cell.node_ids[local]) != contact.node_ids.end();
        }
        for (std::size_t local = 0; local < 3; ++local) {
            const std::size_t next = (local + 1) % 3;
            if (!onContact[local] || !onContact[next])
                continue;
            const std::size_t opposite = (local + 2) % 3;
            const Point2 face0 = meshPoint(mesh, cell.node_ids[local]);
            const Point2 face1 = meshPoint(mesh, cell.node_ids[next]);
            const Point2 oppositePoint =
                meshPoint(mesh, cell.node_ids[opposite]);
            const Point2 tangent = face1 - face0;
            const Real faceLength = tangent.norm();
            if (faceLength <= 1.0e-30)
                return false;
            geometry.faceLocalNodes = {local, next};
            geometry.oppositeLocalNode = opposite;
            geometry.unitNormal =
                Point2{-tangent.y() / faceLength, tangent.x() / faceLength};
            geometry.faceMidpoint = 0.5 * (face0 + face1);
            const Point2 oppositeOffset =
                oppositePoint - geometry.faceMidpoint;
            geometry.oppositeToFaceMidpointDistance = oppositeOffset.norm();
            geometry.normalHeight =
                std::abs(oppositeOffset.dot(geometry.unitNormal));
            if (geometry.oppositeToFaceMidpointDistance <= 1.0e-30 ||
                geometry.normalHeight <= 1.0e-30) {
                return false;
            }
            const Point2 centroid = (face0 + face1 + oppositePoint) / 3.0;
            const Real centroidDistance =
                std::abs((centroid - geometry.faceMidpoint)
                             .dot(geometry.unitNormal));
            geometry.centroidElectricWeight = std::clamp(
                1.0 - centroidDistance / geometry.normalHeight, 0.0, 1.0);
            return true;
        }
    }
    return false;
}

inline bool cellUsesContactElectricFieldFallback(
    const ImpactIonizationModelConfig& config,
    const DeviceMesh& mesh,
    const Cell& cell)
{
    return config.contactElectricFieldFallback &&
           config.drivingForce == "quasi_fermi_gradient" &&
           (config.contactElectricFieldFallbackScope == "contact_boundary_face"
                ? cellTouchesContactBoundaryFace(mesh, cell)
                : cellTouchesContact(mesh, cell));
}

inline Real contactElectricFallbackImpactField(
    const ImpactIonizationModelConfig& config,
    const DeviceMesh& mesh,
    const Cell& cell,
    const Point2& electricGradient,
    Real quasiFermiField,
    const std::function<Real(Index)>& potential,
    Real fieldFactor)
{
    const Real electricMagnitude = electricGradient.norm() * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "cell_gradient_magnitude")
        return electricMagnitude;

    ContactBoundaryFaceGeometry geometry;
    if (!contactBoundaryFaceGeometry(mesh, cell, geometry))
        return quasiFermiField;
    const Real faceNormalField =
        std::abs(electricGradient.dot(geometry.unitNormal)) * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "face_normal")
        return faceNormalField;

    const Index face0 = cell.node_ids[geometry.faceLocalNodes[0]];
    const Index face1 = cell.node_ids[geometry.faceLocalNodes[1]];
    const Index opposite = cell.node_ids[geometry.oppositeLocalNode];
    const Real oneSidedField = std::abs(
        potential(opposite) - 0.5 * (potential(face0) + potential(face1))) /
        geometry.oppositeToFaceMidpointDistance * fieldFactor;
    if (config.contactElectricFieldFallbackMode == "one_sided")
        return oneSidedField;

    return geometry.centroidElectricWeight * faceNormalField +
           (1.0 - geometry.centroidElectricWeight) * quasiFermiField;
}


inline bool edgeTouchesContactElement(const DeviceMesh& mesh,
                                      const std::vector<std::vector<Index>>& edgeCells,
                                      Index edgeId,
                                      const std::vector<bool>& contactNodes)
{
    const Edge& edge = mesh.getEdge(edgeId);
    if (contactNodes[static_cast<std::size_t>(edge.n0)] ||
        contactNodes[static_cast<std::size_t>(edge.n1)]) {
        return true;
    }
    if (edgeId >= edgeCells.size())
        return false;
    for (Index cellId : edgeCells[edgeId]) {
        const Cell& cell = mesh.getCell(cellId);
        for (Index nodeId : cell.node_ids) {
            if (nodeId < mesh.numNodes() &&
                contactNodes[static_cast<std::size_t>(nodeId)]) {
                return true;
            }
        }
    }
    return false;
}

inline Real edgeHighFieldDrivingField(bool qfDrivingForce,
                                      Real qfField,
                                      Real electricField,
                                      const std::vector<std::vector<Index>>& edgeCells,
                                      const DeviceMesh& mesh,
                                      Index edgeId,
                                      const std::vector<bool>& contactNodes)
{
    if (!qfDrivingForce)
        return electricField;
    if (edgeTouchesContactElement(mesh, edgeCells, edgeId, contactNodes))
        return electricField;
    return qfField;
}

inline Real electronQfForAvalancheGradient(Real psi,
                                           Real phin,
                                           Real electronDensity,
                                           Real intrinsicDensity,
                                           Real Vt,
                                           const ImpactIonizationModelConfig& config)
{
    if (!usesQuasiFermiCarrierTruncation(config) || intrinsicDensity <= 0.0)
        return phin;
    const Real carrier = std::max(
        std::max(electronDensity, 0.0),
        config.quasiFermiCarrierTruncation * intrinsicDensity);
    return psi - Vt * std::log(carrier / intrinsicDensity);
}

inline Real holeQfForAvalancheGradient(Real psi,
                                       Real phip,
                                       Real holeDensity,
                                       Real intrinsicDensity,
                                       Real Vt,
                                       const ImpactIonizationModelConfig& config)
{
    if (!usesQuasiFermiCarrierTruncation(config) || intrinsicDensity <= 0.0)
        return phip;
    const Real carrier = std::max(
        std::max(holeDensity, 0.0),
        config.quasiFermiCarrierTruncation * intrinsicDensity);
    return psi + Vt * std::log(carrier / intrinsicDensity);
}

struct EdgeAvalancheDirectionalWeights {
    Real electronNode0 = 0.5;
    Real electronNode1 = 0.5;
    Real holeNode0 = 0.5;
    Real holeNode1 = 0.5;
};

struct CellScalarGradientCache {
    std::vector<Point2> gradients;
    std::vector<Real> areas;
    std::vector<bool> valid;
};

template <typename ValueAt>
inline Point2 cellScalarGradient(
    const DeviceMesh& mesh,
    const Cell&       cell,
    ValueAt&&         valueAt,
    bool&             valid,
    Real&             area)
{
    valid = false;
    area = 0.0;
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
        return Point2::Zero();

    const Index n0 = cell.node_ids[0];
    const Index n1 = cell.node_ids[1];
    const Index n2 = cell.node_ids[2];
    const Node& p0 = mesh.getNode(n0);
    const Node& p1 = mesh.getNode(n1);
    const Node& p2 = mesh.getNode(n2);
    const Real dx10 = p1.x - p0.x;
    const Real dy10 = p1.y - p0.y;
    const Real dx20 = p2.x - p0.x;
    const Real dy20 = p2.y - p0.y;
    const Real det = dx10 * dy20 - dy10 * dx20;
    if (std::abs(det) <= 1.0e-300)
        return Point2::Zero();

    const Real dv10 = valueAt(n1) - valueAt(n0);
    const Real dv20 = valueAt(n2) - valueAt(n0);
    valid = true;
    area = 0.5 * std::abs(det);
    return Point2{
        (dv10 * dy20 - dv20 * dy10) / det,
        (dx10 * dv20 - dx20 * dv10) / det,
    };
}

template <typename ValueAt>
inline CellScalarGradientCache computeCellScalarGradientCache(
    const DeviceMesh& mesh,
    ValueAt&&         valueAt)
{
    CellScalarGradientCache cache;
    cache.gradients.assign(static_cast<std::size_t>(mesh.numCells()), Point2::Zero());
    cache.areas.assign(static_cast<std::size_t>(mesh.numCells()), 0.0);
    cache.valid.assign(static_cast<std::size_t>(mesh.numCells()), false);

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        bool valid = false;
        Real area = 0.0;
        cache.gradients[cellId] = cellScalarGradient(
            mesh, mesh.getCell(cellId), valueAt, valid, area);
        cache.areas[cellId] = area;
        cache.valid[cellId] = valid;
    }
    return cache;
}

template <typename ValueAt>
inline Point2 edgeAveragedCellScalarGradient(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ValueAt&&                              valueAt,
    bool&                                  valid)
{
    valid = false;
    if (edgeId >= edgeCells.size())
        return Point2::Zero();

    Point2 weightedGradient = Point2::Zero();
    Real totalArea = 0.0;
    for (const Index cellId : edgeCells[edgeId]) {
        bool cellValid = false;
        Real area = 0.0;
        const Point2 gradient = cellScalarGradient(
            mesh, mesh.getCell(cellId), valueAt, cellValid, area);
        if (!cellValid || area <= 0.0)
            continue;
        weightedGradient += area * gradient;
        totalArea += area;
    }

    if (totalArea <= 0.0)
        return Point2::Zero();
    valid = true;
    return weightedGradient / totalArea;
}

inline Point2 edgeAveragedCellScalarGradient(
    const std::vector<std::vector<Index>>& edgeCells,
    Index                                  edgeId,
    const CellScalarGradientCache&         cache,
    bool&                                  valid)
{
    valid = false;
    if (edgeId >= edgeCells.size())
        return Point2::Zero();

    Point2 weightedGradient = Point2::Zero();
    Real totalArea = 0.0;
    for (const Index cellId : edgeCells[edgeId]) {
        if (cellId >= cache.valid.size() || !cache.valid[cellId])
            continue;
        const Real area = cache.areas[cellId];
        if (area <= 0.0)
            continue;
        weightedGradient += area * cache.gradients[cellId];
        totalArea += area;
    }

    if (totalArea <= 0.0)
        return Point2::Zero();
    valid = true;
    return weightedGradient / totalArea;
}

inline std::vector<Real> computeNodeCellGradientMagnitudes(
    const std::vector<std::vector<Index>>& nodeCells,
    const CellScalarGradientCache&         cache)
{
    std::vector<Point2> gradients(nodeCells.size(), Point2::Zero());
    for (std::size_t node = 0; node < nodeCells.size(); ++node) {
        Point2 weightedGradient = Point2::Zero();
        Real totalArea = 0.0;
        for (const Index cellId : nodeCells[node]) {
            if (cellId >= cache.valid.size() || !cache.valid[cellId])
                continue;
            const Real area = cache.areas[cellId];
            if (area <= 0.0)
                continue;
            weightedGradient += area * cache.gradients[cellId];
            totalArea += area;
        }
        if (totalArea > 0.0)
            gradients[node] = weightedGradient / totalArea;
    }
    std::vector<Real> fields(nodeCells.size(), 0.0);
    for (std::size_t node = 0; node < nodeCells.size(); ++node)
        fields[node] = gradients[node].norm();
    return fields;
}

template <typename ValueAt>
inline std::vector<Real> computeNodeCellGradientMagnitudes(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    return computeNodeCellGradientMagnitudes(
        nodeCells, computeCellScalarGradientCache(mesh, valueAt));
}

template <typename ValueAt>
inline std::vector<Point2> computeNodeWeightedLeastSquaresGradients(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    std::vector<std::unordered_set<Index>> nodeNeighbors(mesh.numNodes());
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        nodeNeighbors[edge.n0].insert(edge.n1);
        nodeNeighbors[edge.n1].insert(edge.n0);
    }

    std::vector<Point2> fields(mesh.numNodes(), Point2::Zero());
    const CellScalarGradientCache fallbackCache = computeCellScalarGradientCache(
        mesh, [&](Index node) { return valueAt(node); });

    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        const Node& center = mesh.getNode(nodeId);
        const Real centerValue = valueAt(nodeId);
        Real sxx = 0.0;
        Real sxy = 0.0;
        Real syy = 0.0;
        Real sxv = 0.0;
        Real syv = 0.0;

        for (const Index neighborId : nodeNeighbors[nodeId]) {
            const Node& neighbor = mesh.getNode(neighborId);
            const Real dx = neighbor.x - center.x;
            const Real dy = neighbor.y - center.y;
            const Real distance = std::hypot(dx, dy);
            if (distance <= 1.0e-30)
                continue;
            const Real weight = 1.0 / distance;
            const Real dv = valueAt(neighborId) - centerValue;
            sxx += weight * dx * dx;
            sxy += weight * dx * dy;
            syy += weight * dy * dy;
            sxv += weight * dx * dv;
            syv += weight * dy * dv;
        }

        const Real det = sxx * syy - sxy * sxy;
        if (std::abs(det) <= 1.0e-60) {
            Point2 weightedGradient = Point2::Zero();
            Real totalArea = 0.0;
            for (const Index cellId : nodeCells[nodeId]) {
                if (cellId >= fallbackCache.valid.size() || !fallbackCache.valid[cellId])
                    continue;
                const Real area = fallbackCache.areas[cellId];
                if (area <= 0.0)
                    continue;
                weightedGradient += area * fallbackCache.gradients[cellId];
                totalArea += area;
            }
            if (totalArea > 0.0)
                fields[nodeId] = weightedGradient / totalArea;
            continue;
        }

        const Real gradX = (sxv * syy - syv * sxy) / det;
        const Real gradY = (sxx * syv - sxy * sxv) / det;
        fields[nodeId] = Point2{gradX, gradY};
    }
    return fields;
}

template <typename ValueAt>
inline std::vector<Real> computeNodeWeightedLeastSquaresGradientMagnitudes(
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    ValueAt&&                              valueAt)
{
    const std::vector<Point2> gradients = computeNodeWeightedLeastSquaresGradients(
        mesh, nodeCells, std::forward<ValueAt>(valueAt));
    std::vector<Real> fields(gradients.size(), 0.0);
    for (std::size_t node = 0; node < gradients.size(); ++node)
        fields[node] = gradients[node].norm();
    return fields;
}

inline std::vector<Real> computeNodeElectricFields(const VectorXd& psi, const DeviceMesh& mesh, Real fieldFactor)
{
    const std::vector<std::vector<Index>> nodeCells = buildNodeCellMap(mesh);
    std::vector<Real> fields = computeNodeWeightedLeastSquaresGradientMagnitudes(
        mesh, nodeCells, [&](Index node) { return psi(static_cast<int>(node)); });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
}

inline Real edgeQuasiFermiCoefficientField(
    const ImpactIonizationModelConfig&     config,
    Real                                   edgeQfField,
    Real                                   electricField,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const std::vector<bool>&               contactNodes,
    const CellScalarGradientCache&         qfGradientCache,
    Real                                   fieldFactor)
{
    if (usesCellGradientQuasiFermiAvalancheDrive(config)) {
        bool validGradient = false;
        const Point2 gradient = edgeAveragedCellScalarGradient(
            edgeCells, edgeId, qfGradientCache, validGradient);
        return validGradient ? gradient.norm() * fieldFactor : edgeQfField;
    }
    if (config.debugRawVanOverstraeten)
        return edgeQfField;
    return edgeHighFieldDrivingField(
        true, edgeQfField, electricField, edgeCells, mesh, edgeId, contactNodes);
}

inline std::vector<Real> computeElectronAvalancheNodeQuasiFermiDrivingFields(
    const ImpactIonizationModelConfig&     config,
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    const VectorXd&                        psi,
    const VectorXd&                        phin,
    const VectorXd&                        n,
    const std::vector<Real>&               ni,
    Real                                   Vt,
    Real                                   fieldFactor = 1.0)
{
    std::vector<Real> fields = !usesCellGradientQuasiFermiAvalancheDrive(config)
        ? computeNodeScalarGradientMagnitudes(phin, mesh)
        : computeNodeCellGradientMagnitudes(
            mesh, nodeCells, [&](Index node) {
                const int idx = static_cast<int>(node);
                return electronQfForAvalancheGradient(
                    psi(idx), phin(idx), n(idx), ni[node], Vt, config);
            });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
}

inline std::vector<Real> computeHoleAvalancheNodeQuasiFermiDrivingFields(
    const ImpactIonizationModelConfig&     config,
    const DeviceMesh&                      mesh,
    const std::vector<std::vector<Index>>& nodeCells,
    const VectorXd&                        psi,
    const VectorXd&                        phip,
    const VectorXd&                        p,
    const std::vector<Real>&               ni,
    Real                                   Vt,
    Real                                   fieldFactor = 1.0)
{
    std::vector<Real> fields = !usesCellGradientQuasiFermiAvalancheDrive(config)
        ? computeNodeScalarGradientMagnitudes(phip, mesh)
        : computeNodeCellGradientMagnitudes(
            mesh, nodeCells, [&](Index node) {
                const int idx = static_cast<int>(node);
                return holeQfForAvalancheGradient(
                    psi(idx), phip(idx), p(idx), ni[node], Vt, config);
            });
    for (Real& field : fields)
        field *= fieldFactor;
    return fields;
}

template <typename ValueAt>
inline Real edgeMinusGradientUnitDot(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ValueAt&&                              valueAt)
{
    const Edge& edge = mesh.getEdge(edgeId);
    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    if (edge.length <= 1.0e-30)
        return 0.0;

    bool validGradient = false;
    const Point2 gradient = edgeAveragedCellScalarGradient(
        edgeCells, mesh, edgeId, valueAt, validGradient);
    if (!validGradient)
        return 0.0;

    const Point2 minusGradient = -gradient;
    const Real gradientNorm = minusGradient.norm();
    if (gradientNorm <= 1.0e-300)
        return 0.0;

    const Point2 edgeUnit{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
    return std::clamp(edgeUnit.dot(minusGradient / gradientNorm), -1.0, 1.0);
}

inline Real edgeMinusGradientUnitDot(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const CellScalarGradientCache&         cache)
{
    const Edge& edge = mesh.getEdge(edgeId);
    const Node& n0 = mesh.getNode(edge.n0);
    const Node& n1 = mesh.getNode(edge.n1);
    if (edge.length <= 1.0e-30)
        return 0.0;

    bool validGradient = false;
    const Point2 gradient = edgeAveragedCellScalarGradient(
        edgeCells, edgeId, cache, validGradient);
    if (!validGradient)
        return 0.0;

    const Point2 minusGradient = -gradient;
    const Real gradientNorm = minusGradient.norm();
    if (gradientNorm <= 1.0e-300)
        return 0.0;

    const Point2 edgeUnit{(n1.x - n0.x) / edge.length, (n1.y - n0.y) / edge.length};
    return std::clamp(edgeUnit.dot(minusGradient / gradientNorm), -1.0, 1.0);
}

inline EdgeAvalancheDirectionalWeights edgeAvalancheDirectionalWeights(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    const CellScalarGradientCache&         electronGradientCache,
    const CellScalarGradientCache&         holeGradientCache)
{
    EdgeAvalancheDirectionalWeights weights;
    const Real electronDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, electronGradientCache);
    const Real holeDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, holeGradientCache);

    weights.electronNode0 = std::clamp(0.5 + 0.5 * electronDot, 0.0, 1.0);
    weights.electronNode1 = 1.0 - weights.electronNode0;
    weights.holeNode1 = std::clamp(0.5 + 0.5 * holeDot, 0.0, 1.0);
    weights.holeNode0 = 1.0 - weights.holeNode1;
    return weights;
}

template <typename ElectronQfAt, typename HoleQfAt>
inline EdgeAvalancheDirectionalWeights edgeAvalancheDirectionalWeights(
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                      mesh,
    Index                                  edgeId,
    ElectronQfAt&&                         electronQfAt,
    HoleQfAt&&                             holeQfAt)
{
    EdgeAvalancheDirectionalWeights weights;
    const Real electronDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, electronQfAt);
    const Real holeDot = edgeMinusGradientUnitDot(
        edgeCells, mesh, edgeId, holeQfAt);

    weights.electronNode0 = std::clamp(0.5 + 0.5 * electronDot, 0.0, 1.0);
    weights.electronNode1 = 1.0 - weights.electronNode0;
    weights.holeNode1 = std::clamp(0.5 + 0.5 * holeDot, 0.0, 1.0);
    weights.holeNode0 = 1.0 - weights.holeNode1;
    return weights;
}

struct TriangleGssAvalancheSourceRecord {
    Index cellId = 0;
    int localEdge = -1;
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real edgeLength = 0.0;
    Real truncatedPartialVolume = 0.0;
    Real electronCellQfField = 0.0;
    Real holeCellQfField = 0.0;
    Real electronEdgeQfField = 0.0;
    Real holeEdgeQfField = 0.0;
    Real electronMidpointDensity = 0.0;
    Real holeMidpointDensity = 0.0;
    Real electronMobility = 0.0;
    Real holeMobility = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electronFluxProxy = 0.0;
    Real holeFluxProxy = 0.0;
    Real electronSourceIntegral = 0.0;
    Real holeSourceIntegral = 0.0;
    Real edgeSourceIntegral = 0.0;
    Real node0SourceIntegral = 0.0;
    Real node1SourceIntegral = 0.0;
};

inline Real triangleGssEndpointAveragedMobility(
    const MobilityModelConfig&    mobilityConfig,
    const MobilityModel&          mobility,
    const DeviceMesh&             mesh,
    const DopingModel&            doping,
    const std::vector<Material>&  cellMaterials,
    const VectorXd&               n,
    const VectorXd&               p,
    Index                         cellId,
    Index                         node0,
    Index                         node1,
    CarrierType                   carrier,
    Real                          drivingField)
{
    const Material& material = cellMaterials.at(static_cast<std::size_t>(cellId));
    const auto atNode = [&](Index node) {
        const Real mobilityDoping = nodeMobilityDopingConcentration(
            mesh, doping, node, cellId, &mobilityConfig);
        return carrier == CarrierType::Electron
            ? mobility.electronMobility(
                material, mobilityDoping, n(static_cast<int>(node)),
                p(static_cast<int>(node)), drivingField, 0.0)
            : mobility.holeMobility(
                material, mobilityDoping, n(static_cast<int>(node)),
                p(static_cast<int>(node)), drivingField, 0.0);
    };
    return 0.5 * (atNode(node0) + atNode(node1));
}

inline std::vector<TriangleGssAvalancheSourceRecord>
triangleGssAvalancheSourceRecordsForCell(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<Index>&          cellEdgeIds,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              cellId,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "triangle GSS avalanche source requires Tri3 cells; unsupported cell " +
            std::to_string(cellId));
    }

    bool electronGradientValid = false;
    bool holeGradientValid = false;
    Real electronDoubleArea = 0.0;
    Real holeDoubleArea = 0.0;
    const Point2 electronGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phin(static_cast<int>(node));
        }, electronGradientValid, electronDoubleArea);
    const Point2 holeGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phip(static_cast<int>(node));
        }, holeGradientValid, holeDoubleArea);
    if (!electronGradientValid || !holeGradientValid)
        return {};

    const Real electronCellField = electronGradient.norm() * fieldFactor;
    const Real holeCellField = holeGradient.norm() * fieldFactor;
    Real electronImpactField = electronCellField;
    Real holeImpactField = holeCellField;
    if (cellUsesContactElectricFieldFallback(config, mesh, cell)) {
        bool electricGradientValid = false;
        Real electricDoubleArea = 0.0;
        const Point2 electricGradient =
            cellScalarGradient(mesh, cell, [&](Index node) {
                return psi(static_cast<int>(node));
            }, electricGradientValid, electricDoubleArea);
        if (!electricGradientValid)
            return {};
        const auto potential = [&](Index node) {
            return psi(static_cast<int>(node));
        };
        electronImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, electronCellField,
            potential, fieldFactor);
        holeImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, holeCellField,
            potential, fieldFactor);
    }
    const Real electronAlpha = impact.electronCoefficient(electronImpactField);
    const Real holeAlpha = impact.holeCoefficient(holeImpactField);
    std::vector<TriangleGssAvalancheSourceRecord> records;
    records.reserve(3);

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const Index node0 = cell.node_ids[static_cast<std::size_t>(localEdge)];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        Index edgeId = mesh.numEdges();
        for (Index candidate : cellEdgeIds) {
            const Edge& edge = mesh.getEdge(candidate);
            if ((edge.n0 == node0 && edge.n1 == node1) ||
                (edge.n0 == node1 && edge.n1 == node0)) {
                edgeId = candidate;
                break;
            }
        }
        if (edgeId >= mesh.numEdges())
            throw std::runtime_error("triangle GSS avalanche source could not map a cell edge");
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.length <= 1.0e-30)
            continue;

        TriangleGssAvalancheSourceRecord record;
        record.cellId = cellId;
        record.localEdge = localEdge;
        record.edgeId = edgeId;
        record.node0 = node0;
        record.node1 = node1;
        record.edgeLength = edge.length;
        record.truncatedPartialVolume =
            geniusTri3TruncatedPartialVolumeWithEdge(
                mesh, cell, node0, node1) * config.sourceGeometryScale;
        record.electronCellQfField = electronCellField;
        record.holeCellQfField = holeCellField;
        record.electronEdgeQfField =
            std::abs(phin(static_cast<int>(node1)) - phin(static_cast<int>(node0))) /
            edge.length * fieldFactor;
        record.holeEdgeQfField =
            std::abs(phip(static_cast<int>(node1)) - phip(static_cast<int>(node0))) /
            edge.length * fieldFactor;
        record.electronMidpointDensity = gssElectronAvalancheMidpointDensity(
            n(static_cast<int>(node0)), n(static_cast<int>(node1)),
            psi(static_cast<int>(node0)), psi(static_cast<int>(node1)), Vt);
        record.holeMidpointDensity = gssHoleAvalancheMidpointDensity(
            p(static_cast<int>(node0)), p(static_cast<int>(node1)),
            psi(static_cast<int>(node0)), psi(static_cast<int>(node1)), Vt);
        record.electronMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Electron, record.electronEdgeQfField);
        record.holeMobility = triangleGssEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId, node0, node1,
            CarrierType::Hole, record.holeEdgeQfField);
        record.electronAlpha = electronAlpha;
        record.holeAlpha = holeAlpha;
        record.electronFluxProxy = record.electronMobility *
            record.electronMidpointDensity * record.electronEdgeQfField;
        record.holeFluxProxy = record.holeMobility *
            record.holeMidpointDensity * record.holeEdgeQfField;
        record.electronSourceIntegral = electronAlpha *
            record.electronFluxProxy * record.truncatedPartialVolume;
        record.holeSourceIntegral = holeAlpha *
            record.holeFluxProxy * record.truncatedPartialVolume;
        record.edgeSourceIntegral =
            record.electronSourceIntegral + record.holeSourceIntegral;
        record.node0SourceIntegral = 0.5 * record.edgeSourceIntegral;
        record.node1SourceIntegral = 0.5 * record.edgeSourceIntegral;
        records.push_back(record);
    }
    return records;
}

inline std::vector<TriangleGssAvalancheSourceRecord>
triangleGssAvalancheSourceRecords(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    if (isSurfaceMobilityModel(mobilityConfig)) {
        throw std::invalid_argument(
            "triangle GSS avalanche source does not support surface mobility");
    }
    std::vector<TriangleGssAvalancheSourceRecord> records;
    const auto cellEdges = buildCellEdgeMap(edgeCells, mesh);
    records.reserve(static_cast<std::size_t>(mesh.numCells()) * 3);

    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        auto cellRecords = triangleGssAvalancheSourceRecordsForCell(
            config, impact, mobilityConfig, mobility,
            cellEdges.at(static_cast<std::size_t>(cellId)),
            mesh, doping, cellMaterials, cellId, psi, phin, phip, n, p, Vt,
            fieldFactor);
        records.insert(records.end(), cellRecords.begin(), cellRecords.end());
    }
    return records;
}

#include "vela/equation/ElementEdgeGssLauxAD.inl"

struct ElementEdgeGssLauxAvalancheSourceRecord {
    Index cellId = 0;
    std::array<Index, 3> edgeIds{};
    std::array<Real, 3> edgeLengths{};
    std::array<Real, 3> edgePartialVolumes{};
    std::array<Real, 3> vertexMeasures{};
    std::array<Real, 3> electronMobilities{};
    std::array<Real, 3> holeMobilities{};
    std::array<Real, 3> electronSignedEdgeFlux{};
    std::array<Real, 3> holeSignedEdgeFlux{};
    Point2 electronCurrentVector = Point2::Zero();
    Point2 holeCurrentVector = Point2::Zero();
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    std::array<Real, 3> electronSourceIntegrals{};
    std::array<Real, 3> holeSourceIntegrals{};
    std::array<Real, 3> combinedSourceIntegrals{};
};

inline ElementEdgeGssLauxAvalancheSourceRecord
elementEdgeGssLauxAvalancheSourceRecordForCell(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<Index>&          cellEdgeIds,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              cellId,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    if (!usesElementEdgeGssLauxAvalancheSource(config)) {
        throw std::invalid_argument(
            "element-edge GSS/Laux record requires its canonical configuration");
    }
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge GSS/Laux avalanche source requires Tri3 cells");
    }

    bool electricGradientValid = false;
    Real electricDoubleArea = 0.0;
    const Point2 electricGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return psi(static_cast<int>(node));
        }, electricGradientValid, electricDoubleArea);
    bool electronGradientValid = false;
    bool holeGradientValid = false;
    Real electronDoubleArea = 0.0;
    Real holeDoubleArea = 0.0;
    const Point2 electronGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phin(static_cast<int>(node));
        }, electronGradientValid, electronDoubleArea);
    const Point2 holeGradient =
        cellScalarGradient(mesh, cell, [&](Index node) {
            return phip(static_cast<int>(node));
        }, holeGradientValid, holeDoubleArea);
    if (!electricGradientValid || !electronGradientValid ||
        !holeGradientValid) {
        throw std::invalid_argument(
            "degenerate triangle cannot evaluate avalanche driving fields");
    }

    ElementEdgeGssLauxAvalancheSourceRecord record;
    record.cellId = cellId;
    record.edgePartialVolumes =
        tri3ElementEdgeBoxPartialVolumes(mesh, cell);
    record.vertexMeasures =
        tri3ElementVertexBoxMeasures(mesh, cell);
    const bool useContactFallback =
        cellUsesContactElectricFieldFallback(config, mesh, cell);
    const Real electronQfField = electronGradient.norm() * fieldFactor;
    const Real holeQfField = holeGradient.norm() * fieldFactor;
    if (config.drivingForce == "electric_field") {
        record.electronImpactField = electricGradient.norm() * fieldFactor;
        record.holeImpactField = record.electronImpactField;
    } else if (useContactFallback) {
        const auto potential = [&](Index node) {
            return psi(static_cast<int>(node));
        };
        record.electronImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, electronQfField,
            potential, fieldFactor);
        record.holeImpactField = contactElectricFallbackImpactField(
            config, mesh, cell, electricGradient, holeQfField,
            potential, fieldFactor);
    } else {
        record.electronImpactField = electronQfField;
        record.holeImpactField = holeQfField;
    }
    record.electronAlpha =
        impact.electronCoefficient(record.electronImpactField);
    record.holeAlpha =
        impact.holeCoefficient(record.holeImpactField);

    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const std::size_t local = static_cast<std::size_t>(localEdge);
        const Index node0 = cell.node_ids[local];
        const Index node1 =
            cell.node_ids[static_cast<std::size_t>((localEdge + 1) % 3)];
        const Index edgeId =
            edgeIdForNodePair(mesh, cellEdgeIds, node0, node1);
        if (edgeId >= mesh.numEdges()) {
            throw std::runtime_error(
                "element-edge GSS/Laux source could not map a cell edge");
        }
        const Real edgeLength =
            (meshPoint(mesh, node1) - meshPoint(mesh, node0)).norm();
        if (edgeLength <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot evaluate an SG current");
        }
        const Real electricEdgeField =
            std::abs(
                psi(static_cast<int>(node1)) -
                psi(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const bool qfMobility =
            mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
        record.edgeIds[local] = edgeId;
        record.edgeLengths[local] = edgeLength;

        const Real electronEdgeField =
            std::abs(
                phin(static_cast<int>(node1)) -
                phin(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const Real holeEdgeField =
            std::abs(
                phip(static_cast<int>(node1)) -
                phip(static_cast<int>(node0))) /
            edgeLength * fieldFactor;
        const Real electronMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Electron,
                qfMobility ? electronEdgeField : electricEdgeField);
        const Real holeMobility =
            triangleGssEndpointAveragedMobility(
                mobilityConfig, mobility, mesh, doping, cellMaterials, n, p, cellId,
                node0, node1, CarrierType::Hole,
                qfMobility ? holeEdgeField : electricEdgeField);
        record.electronMobilities[local] = electronMobility;
        record.holeMobilities[local] = holeMobility;
        record.electronSignedEdgeFlux[local] =
            electronMobility > 0.0
            ? sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni.at(node0), ni.at(node1),
                psi(static_cast<int>(node0)), psi(static_cast<int>(node1)),
                phin(static_cast<int>(node0)), phin(static_cast<int>(node1)),
                Vt, electronMobility * Vt * fieldFactor / edgeLength)
            : 0.0;
        record.holeSignedEdgeFlux[local] =
            holeMobility > 0.0
            ? sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni.at(node0), ni.at(node1),
                psi(static_cast<int>(node0)), psi(static_cast<int>(node1)),
                phip(static_cast<int>(node0)), phip(static_cast<int>(node1)),
                Vt, holeMobility * Vt * fieldFactor / edgeLength)
            : 0.0;
    }

    record.electronCurrentVector = gssLauxTri3CurrentVector(
        mesh, cell, record.electronSignedEdgeFlux);
    record.holeCurrentVector = gssLauxTri3CurrentVector(
        mesh, cell, record.holeSignedEdgeFlux);
    std::array<Real, 3> localPsi{};
    std::array<Real, 3> localPhin{};
    std::array<Real, 3> localPhip{};
    std::array<Real, 3> localElectronDensity{};
    std::array<Real, 3> localHoleDensity{};
    std::array<Real, 3> localIntrinsicDensity{};
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        const Index node = cell.node_ids[localNode];
        localPsi[localNode] = psi(static_cast<int>(node));
        localPhin[localNode] = phin(static_cast<int>(node));
        localPhip[localNode] = phip(static_cast<int>(node));
        localElectronDensity[localNode] = n(static_cast<int>(node));
        localHoleDensity[localNode] = p(static_cast<int>(node));
        localIntrinsicDensity[localNode] = ni.at(node);
    }
    const auto sharedSourceIntegrals =
        elementEdgeGssLauxAvalancheSourceIntegralsLocal<Real>(
            config, mobilityConfig, mobility, cellEdgeIds, mesh, doping,
            cellMaterials, cellId, localPsi, localPhin, localPhip,
            localElectronDensity, localHoleDensity, localIntrinsicDensity,
            Vt, fieldFactor);
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        record.electronSourceIntegrals[localNode] =
            sharedSourceIntegrals.electron[localNode];
        record.holeSourceIntegrals[localNode] =
            sharedSourceIntegrals.hole[localNode];
        record.combinedSourceIntegrals[localNode] =
            sharedSourceIntegrals.combined[localNode];
    }
    return record;
}

inline std::vector<ElementEdgeGssLauxAvalancheSourceRecord>
elementEdgeGssLauxAvalancheSourceRecords(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    if (isSurfaceMobilityModel(mobilityConfig)) {
        throw std::invalid_argument(
            "element-edge GSS/Laux avalanche source does not support surface mobility");
    }
    const auto cellEdges = buildCellEdgeMap(edgeCells, mesh);
    std::vector<ElementEdgeGssLauxAvalancheSourceRecord> records;
    records.reserve(static_cast<std::size_t>(mesh.numCells()));
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        records.push_back(
            elementEdgeGssLauxAvalancheSourceRecordForCell(
                config, impact, mobilityConfig, mobility,
                cellEdges.at(static_cast<std::size_t>(cellId)),
                mesh, doping, cellMaterials, cellId, psi, phin, phip, n, p,
                ni, Vt, fieldFactor));
    }
    return records;
}

struct SgEdgeCurrentAvalancheSourceRecord {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Real edgeLength = 0.0;
    Real edgeCouple = 0.0;
    Real edgeAreaProxy = 0.0;
    Real electricField = 0.0;
    Real electronImpactField = 0.0;
    Real holeImpactField = 0.0;
    Real electronAlpha = 0.0;
    Real holeAlpha = 0.0;
    Real electronMobility = 0.0;
    Real holeMobility = 0.0;
    Real electronRawFluxProxy = 0.0;
    Real holeRawFluxProxy = 0.0;
    Real electronRawSignedFluxProxy = 0.0;
    Real holeRawSignedFluxProxy = 0.0;
    Real electronReconstructedFluxProxy = 0.0;
    Real holeReconstructedFluxProxy = 0.0;
    Real electronFluxProxy = 0.0;
    Real holeFluxProxy = 0.0;
    Real electronFinalOverRawFluxProxy = 0.0;
    Real holeFinalOverRawFluxProxy = 0.0;
    Real electronSourceIntegral = 0.0;
    Real holeSourceIntegral = 0.0;
    Real edgeSourceIntegral = 0.0;
    Real electronNode0SourceIntegral = 0.0;
    Real electronNode1SourceIntegral = 0.0;
    Real holeNode0SourceIntegral = 0.0;
    Real holeNode1SourceIntegral = 0.0;
    Real node0SourceIntegral = 0.0;
    Real node1SourceIntegral = 0.0;
    SgElectronVariableNiFluxDecomposition electronSgFluxDecomposition;
    Real electronSgProductionSignedFluxNative = 0.0;
    Real electronSgReconstructionRelativeError = 0.0;
    Real electronSgProductionVsReferenceRelativeError = 0.0;
    Real electronSgStableVsReferenceRelativeError = 0.0;
    bool electronSgDiagnosticsCollected = false;
};

struct SgAvalancheSourceComponentIntegrals {
    std::vector<Real> electron;
    std::vector<Real> hole;
    std::vector<Real> combined;
};

inline std::vector<SgEdgeCurrentAvalancheSourceRecord> sgEdgeCurrentAvalancheSourceRecords(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0,
    bool                               collectElectronSgDiagnostics = false)
{
    std::vector<SgEdgeCurrentAvalancheSourceRecord> records;
    records.reserve(mesh.numEdges());
    const bool qfImpact = usesQuasiFermiAvalancheDrivingForce(config);
    const bool currentAlignedImpact =
        !config.debugRawVanOverstraeten && usesCurrentAlignedAvalancheDrivingForce(config);
    const bool cellCurrentReconstructedCurrent = usesCellCurrentReconstructedAvalancheCurrent(config);
    const bool cellVectorCurrentReconstructedCurrent = usesCellVectorCurrentReconstructedAvalancheCurrent(config);
    const bool dualFaceVectorCurrentMagnitude = config.currentMagnitudeMode == "dual_face_vector_mag";
    const bool usesReconstructedSgCurrent = cellCurrentReconstructedCurrent || cellVectorCurrentReconstructedCurrent || dualFaceVectorCurrentMagnitude;
    const bool directionalEdgePartition = usesDirectionalEdgeAvalancheSourcePartition(config);
    const bool qfMobility = mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
    const std::vector<bool> contactNodes = contactNodeMask(mesh);
    const CellScalarGradientCache electronQfGradientCache = computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const int idx = static_cast<int>(node);
            return electronQfForAvalancheGradient(
                psi(idx), phin(idx), n(idx), ni[node], Vt, config);
        });
    const CellScalarGradientCache holeQfGradientCache = computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const int idx = static_cast<int>(node);
            return holeQfForAvalancheGradient(
                psi(idx), phip(idx), p(idx), ni[node], Vt, config);
        });

    std::vector<Real> rawElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawSignedElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> rawSignedHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> reconstructedElectronFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    std::vector<Real> reconstructedHoleFlux(static_cast<std::size_t>(mesh.numEdges()), 0.0);
    if (usesReconstructedSgCurrent) {
        for (Index e = 0; e < mesh.numEdges(); ++e) {
            const Edge& edge = mesh.getEdge(e);
            const Real h = edge.length;
            if (h <= 1.0e-30 || edge.couple <= 0.0)
                continue;
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            const Real psi_i = psi(i);
            const Real psi_j = psi(j);
            const Real phin_i = phin(i);
            const Real phin_j = phin(j);
            const Real phip_i = phip(i);
            const Real phip_j = phip(j);
            const Real electronQf_i = electronQfForAvalancheGradient(
                psi_i, phin_i, n(i), ni[edge.n0], Vt, config);
            const Real electronQf_j = electronQfForAvalancheGradient(
                psi_j, phin_j, n(j), ni[edge.n1], Vt, config);
            const Real holeQf_i = holeQfForAvalancheGradient(
                psi_i, phip_i, p(i), ni[edge.n0], Vt, config);
            const Real holeQf_j = holeQfForAvalancheGradient(
                psi_j, phip_j, p(j), ni[edge.n1], Vt, config);
            const Real electricField = std::abs((psi_j - psi_i) / h) * fieldFactor;
            const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h) * fieldFactor;
            const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h) * fieldFactor;
            const Real electronMobilityField = qfMobility ? electronQfField : electricField;
            const Real holeMobilityField = qfMobility ? holeQfField : electricField;
            const Real mun = edgeMobility(
                edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Electron,
                electronMobilityField, &mobilityConfig, &psi);
            if (mun > 0.0) {
                const Real signedFlux = sgElectronContinuityFluxFromQuasiFermiVariableNi(
                    ni[edge.n0], ni[edge.n1], psi_i, psi_j, phin_i, phin_j,
                    Vt, mun * Vt * fieldFactor / h);
                rawSignedElectronFlux[static_cast<std::size_t>(e)] = signedFlux;
                rawElectronFlux[static_cast<std::size_t>(e)] = std::abs(signedFlux);
            }
            const Real mup = edgeMobility(
                edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Hole,
                holeMobilityField, &mobilityConfig, &psi);
            if (mup > 0.0) {
                const Real signedFlux = sgHoleContinuityFluxFromQuasiFermiVariableNi(
                    ni[edge.n0], ni[edge.n1], psi_i, psi_j, phip_i, phip_j,
                    Vt, mup * Vt * fieldFactor / h);
                rawSignedHoleFlux[static_cast<std::size_t>(e)] = signedFlux;
                rawHoleFlux[static_cast<std::size_t>(e)] = std::abs(signedFlux);
            }
        }

        const std::vector<std::vector<Index>> cellEdges = buildCellEdgeMap(edgeCells, mesh);
        for (Index e = 0; e < mesh.numEdges(); ++e) {
            reconstructedElectronFlux[static_cast<std::size_t>(e)] = dualFaceVectorCurrentMagnitude
                ? medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                    e, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
                : (cellVectorCurrentReconstructedCurrent
                    ? cellVectorReconstructedEdgeFluxMagnitude(
                        e, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
                    : cellSmoothedEdgeFluxMagnitude(e, rawElectronFlux, edgeCells, cellEdges));
            reconstructedHoleFlux[static_cast<std::size_t>(e)] = dualFaceVectorCurrentMagnitude
                ? medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                    e, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
                : (cellVectorCurrentReconstructedCurrent
                    ? cellVectorReconstructedEdgeFluxMagnitude(
                        e, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
                    : cellSmoothedEdgeFluxMagnitude(e, rawHoleFlux, edgeCells, cellEdges));
        }
    }

    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Real h = edge.length;
        if (h <= 1.0e-30 || edge.couple <= 0.0)
            continue;

        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real psi_i = psi(i);
        const Real psi_j = psi(j);
        const Real phin_i = phin(i);
        const Real phin_j = phin(j);
        const Real phip_i = phip(i);
        const Real phip_j = phip(j);

        const Real electronQf_i = electronQfForAvalancheGradient(
            psi_i, phin_i, n(i), ni[edge.n0], Vt, config);
        const Real electronQf_j = electronQfForAvalancheGradient(
            psi_j, phin_j, n(j), ni[edge.n1], Vt, config);
        const Real holeQf_i = holeQfForAvalancheGradient(
            psi_i, phip_i, p(i), ni[edge.n0], Vt, config);
        const Real holeQf_j = holeQfForAvalancheGradient(
            psi_j, phip_j, p(j), ni[edge.n1], Vt, config);
        const Real electricField = std::abs((psi_j - psi_i) / h) * fieldFactor;
        const Real electronQfField = std::abs((electronQf_j - electronQf_i) / h) * fieldFactor;
        const Real holeQfField = std::abs((holeQf_j - holeQf_i) / h) * fieldFactor;
        const Real electronCoefficientField = qfImpact
            ? edgeQuasiFermiCoefficientField(
                config, electronQfField, electricField, edgeCells, mesh, e,
                contactNodes, electronQfGradientCache, fieldFactor)
            : electricField;
        const Real holeCoefficientField = qfImpact
            ? edgeQuasiFermiCoefficientField(
                config, holeQfField, electricField, edgeCells, mesh, e,
                contactNodes, holeQfGradientCache, fieldFactor)
            : electricField;
        const Real electronMobilityField = qfMobility ? electronQfField : electricField;
        const Real holeMobilityField = qfMobility ? holeQfField : electricField;

        const Real nAvg = 0.5 * (n(i) + n(j));
        const Real pAvg = 0.5 * (p(i) + p(j));
        const Real nMid = cellReconstructedAvalancheMidpointDensity(
            config, n(i), n(j), psi_i, psi_j, Vt);
        const Real pMid = cellReconstructedAvalancheMidpointDensity(
            config, p(i), p(j), psi_j, psi_i, Vt);
        const Real signedElectricField01 = -(psi_j - psi_i) / h * fieldFactor;

        const Real edgeArea = avalancheSourceEdgeArea(config, edgeCells, mesh, e);
        SgEdgeCurrentAvalancheSourceRecord record;
        record.edgeId = e;
        record.node0 = edge.n0;
        record.node1 = edge.n1;
        record.edgeLength = h;
        record.edgeCouple = edge.couple;
        record.edgeAreaProxy = edgeArea;
        record.electricField = electricField;

        const Real mun = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Electron,
            electronMobilityField, &mobilityConfig, &psi);
        record.electronMobility = mun;
        const Real mup = edgeMobility(
            edgeCells, mesh, doping, mobility, cellMaterials, e, CarrierType::Hole,
            holeMobilityField, &mobilityConfig, &psi);
        record.holeMobility = mup;
        constexpr bool IncludeElectronNiGradientDrift = true;
        const Real electronContinuityCoefficient =
            mun > 0.0 ? mun * Vt * fieldFactor / h : 0.0;
        const Real electronContinuityFlux01 = mun > 0.0
            ? sgElectronContinuityFluxFromQuasiFermiVariableNi(
                ni[edge.n0],
                ni[edge.n1],
                psi_i,
                psi_j,
                phin_i,
                phin_j,
                Vt,
                electronContinuityCoefficient,
                IncludeElectronNiGradientDrift)
            : 0.0;
        const Real holeContinuityFlux01 = mup > 0.0
            ? sgHoleContinuityFluxFromQuasiFermiVariableNi(
                ni[edge.n0],
                ni[edge.n1],
                psi_i,
                psi_j,
                phip_i,
                phip_j,
                Vt,
                mup * Vt * fieldFactor / h)
            : 0.0;
        // SG continuity fluxes use the particle-flux convention, so the
        // physical charge-current magnitude is |F_p - F_n|.
        const Real conservedTotalFluxMagnitude =
            conservedTotalCurrentFluxMagnitude(
                electronContinuityFlux01, holeContinuityFlux01);
        if (collectElectronSgDiagnostics) {
            record.electronSgDiagnosticsCollected = true;
            record.electronSgFluxDecomposition =
                sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
                    ni[edge.n0],
                    ni[edge.n1],
                    psi_i,
                    psi_j,
                    phin_i,
                    phin_j,
                    Vt,
                    electronContinuityCoefficient,
                    IncludeElectronNiGradientDrift);
            record.electronSgProductionSignedFluxNative = electronContinuityFlux01;
            const Real electronSgReconstructionScale = std::max({
                std::abs(record.electronSgProductionSignedFluxNative),
                std::abs(record.electronSgFluxDecomposition.reconstructedFlux),
                std::abs(record.electronSgFluxDecomposition.coef)
                    * (std::abs(record.electronSgFluxDecomposition.leftTerm)
                       + std::abs(record.electronSgFluxDecomposition.rightTerm)),
                Real{1.0e-300}});
            record.electronSgReconstructionRelativeError =
                std::abs(record.electronSgProductionSignedFluxNative
                         - record.electronSgFluxDecomposition.reconstructedFlux)
                / electronSgReconstructionScale;

            const Real referenceScale = std::max({
                std::abs(record.electronSgProductionSignedFluxNative),
                record.electronSgFluxDecomposition.highPrecisionReferenceTermScale,
                Real{1.0e-300}});
            const long double referenceScaleLong =
                static_cast<long double>(referenceScale);
            const long double referenceFluxLong = static_cast<long double>(
                record.electronSgFluxDecomposition.highPrecisionReferenceFlux);
            record.electronSgProductionVsReferenceRelativeError =
                static_cast<Real>(std::abs(
                    static_cast<long double>(
                        record.electronSgProductionSignedFluxNative)
                    - referenceFluxLong) / referenceScaleLong);
            record.electronSgStableVsReferenceRelativeError =
                static_cast<Real>(std::abs(
                    static_cast<long double>(
                        record.electronSgFluxDecomposition.stableFactorizedFlux)
                    - referenceFluxLong) / referenceScaleLong);
        }

        if (mun > 0.0) {
            record.electronImpactField = currentAlignedImpact
                ? parallelCurrentAvalancheDrivingField(
                    signedElectricField01, electronContinuityFlux01)
                : electronAvalancheDrivingField(
                    config, electronCoefficientField, electricField, nAvg);
            record.electronRawSignedFluxProxy = electronContinuityFlux01;
            record.electronRawFluxProxy = std::abs(electronContinuityFlux01);
            record.electronReconstructedFluxProxy = usesReconstructedSgCurrent
                ? reconstructedElectronFlux[static_cast<std::size_t>(e)]
                : record.electronRawFluxProxy;
            record.electronFluxProxy = selectAvalancheCurrentFluxProxy(
                config,
                record.electronRawFluxProxy,
                record.electronReconstructedFluxProxy,
                mun,
                nMid,
                record.electronImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            record.electronFinalOverRawFluxProxy = record.electronRawFluxProxy > 0.0
                ? record.electronFluxProxy / record.electronRawFluxProxy
                : 0.0;
            record.electronAlpha = impact.electronCoefficient(record.electronImpactField);
            record.electronSourceIntegral =
                record.electronAlpha * record.electronFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.electronSourceIntegral;
        }

        if (mup > 0.0) {
            record.holeImpactField = currentAlignedImpact
                ? parallelCurrentAvalancheDrivingField(
                    signedElectricField01, holeContinuityFlux01)
                : holeAvalancheDrivingField(
                    config, holeCoefficientField, electricField, pAvg);
            record.holeRawSignedFluxProxy = holeContinuityFlux01;
            record.holeRawFluxProxy = std::abs(holeContinuityFlux01);
            record.holeReconstructedFluxProxy = usesReconstructedSgCurrent
                ? reconstructedHoleFlux[static_cast<std::size_t>(e)]
                : record.holeRawFluxProxy;
            record.holeFluxProxy = selectAvalancheCurrentFluxProxy(
                config,
                record.holeRawFluxProxy,
                record.holeReconstructedFluxProxy,
                mup,
                pMid,
                record.holeImpactField,
                electricField,
                conservedTotalFluxMagnitude);
            record.holeFinalOverRawFluxProxy = record.holeRawFluxProxy > 0.0
                ? record.holeFluxProxy / record.holeRawFluxProxy
                : 0.0;
            record.holeAlpha = impact.holeCoefficient(record.holeImpactField);
            record.holeSourceIntegral =
                record.holeAlpha * record.holeFluxProxy * edgeArea;
            record.edgeSourceIntegral += record.holeSourceIntegral;
        }

        EdgeAvalancheDirectionalWeights weights;
        if (directionalEdgePartition) {
            weights = edgeAvalancheDirectionalWeights(
                edgeCells,
                mesh,
                e,
                electronQfGradientCache,
                holeQfGradientCache);
        }
        record.electronNode0SourceIntegral =
            weights.electronNode0 * record.electronSourceIntegral;
        record.electronNode1SourceIntegral =
            weights.electronNode1 * record.electronSourceIntegral;
        record.holeNode0SourceIntegral =
            weights.holeNode0 * record.holeSourceIntegral;
        record.holeNode1SourceIntegral =
            weights.holeNode1 * record.holeSourceIntegral;
        record.node0SourceIntegral =
            record.electronNode0SourceIntegral + record.holeNode0SourceIntegral;
        record.node1SourceIntegral =
            record.electronNode1SourceIntegral + record.holeNode1SourceIntegral;
        records.push_back(record);
    }
    return records;
}

inline void addCellMappedEdgeSourceToNodes(
    std::vector<Real>&                       target,
    const std::vector<std::vector<Index>>&   edgeCells,
    const DeviceMesh&                        mesh,
    const SgEdgeCurrentAvalancheSourceRecord& record,
    Real                                     sourceIntegral)
{
    if (record.edgeId >= edgeCells.size() || sourceIntegral == 0.0)
        return;
    const auto& cells = edgeCells[record.edgeId];
    Real areaSum = 0.0;
    for (Index cellId : cells) {
        if (cellId < mesh.numCells())
            areaSum += triangleArea(mesh, mesh.getCell(cellId));
    }
    if (areaSum <= 0.0) {
        if (record.node0 < target.size()) target[record.node0] += 0.5 * sourceIntegral;
        if (record.node1 < target.size()) target[record.node1] += 0.5 * sourceIntegral;
        return;
    }
    for (Index cellId : cells) {
        if (cellId >= mesh.numCells())
            continue;
        const Cell& cell = mesh.getCell(cellId);
        if (cell.node_ids.empty())
            continue;
        const Real cellShare = sourceIntegral * triangleArea(mesh, cell) / areaSum;
        const Real nodeShare = cellShare / static_cast<Real>(cell.node_ids.size());
        for (Index nodeId : cell.node_ids) {
            if (nodeId < target.size())
                target[nodeId] += nodeShare;
        }
    }
}

inline void addMappedEdgeSourceToNodes(
    const ImpactIonizationModelConfig&       config,
    std::vector<Real>&                       target,
    const std::vector<std::vector<Index>>&   edgeCells,
    const DeviceMesh&                        mesh,
    const SgEdgeCurrentAvalancheSourceRecord& record,
    Real                                     node0DirectionalSource,
    Real                                     node1DirectionalSource,
    Real                                     sourceIntegral)
{
    if (config.sourceMappingMode == "cell_F_cell_alpha_cell_G_to_node") {
        addCellMappedEdgeSourceToNodes(target, edgeCells, mesh, record, sourceIntegral);
        return;
    }
    if (config.sourceMappingMode == "edge_F_edge_alpha_edge_G_to_node") {
        if (record.node0 < target.size()) target[record.node0] += 0.5 * sourceIntegral;
        if (record.node1 < target.size()) target[record.node1] += 0.5 * sourceIntegral;
        return;
    }
    if (record.node0 < target.size()) target[record.node0] += node0DirectionalSource;
    if (record.node1 < target.size()) target[record.node1] += node1DirectionalSource;
}
inline SgAvalancheSourceComponentIntegrals sgEdgeCurrentAvalancheSourceComponentIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    SgAvalancheSourceComponentIntegrals source;
    source.electron.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.hole.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.combined.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    const auto records = sgEdgeCurrentAvalancheSourceRecords(
        config,
        impact,
        mobilityConfig,
        mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt,
        fieldFactor);
    for (const auto& record : records) {
        addMappedEdgeSourceToNodes(
            config, source.electron, edgeCells, mesh, record,
            record.electronNode0SourceIntegral,
            record.electronNode1SourceIntegral,
            record.electronSourceIntegral);
        addMappedEdgeSourceToNodes(
            config, source.hole, edgeCells, mesh, record,
            record.holeNode0SourceIntegral,
            record.holeNode1SourceIntegral,
            record.holeSourceIntegral);
        addMappedEdgeSourceToNodes(
            config, source.combined, edgeCells, mesh, record,
            record.node0SourceIntegral,
            record.node1SourceIntegral,
            record.edgeSourceIntegral);
    }
    return source;
}

inline SgAvalancheSourceComponentIntegrals
currentDensityAvalancheSourceComponentIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    if (usesElementEdgeGssLauxAvalancheSource(config)) {
        SgAvalancheSourceComponentIntegrals source;
        source.electron.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        source.hole.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        source.combined.assign(
            static_cast<std::size_t>(mesh.numNodes()), 0.0);
        const auto records = elementEdgeGssLauxAvalancheSourceRecords(
            config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
            cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor);
        for (const auto& record : records) {
            const Cell& cell = mesh.getCell(record.cellId);
            for (std::size_t localNode = 0; localNode < 3; ++localNode) {
                const Index node = cell.node_ids[localNode];
                source.electron[node] +=
                    record.electronSourceIntegrals[localNode];
                source.hole[node] +=
                    record.holeSourceIntegrals[localNode];
                source.combined[node] +=
                    record.combinedSourceIntegrals[localNode];
            }
        }
        return source;
    }

    if (!usesTriangleGssAvalancheSource(config)) {
        return sgEdgeCurrentAvalancheSourceComponentIntegrals(
            config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
            cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor);
    }

    SgAvalancheSourceComponentIntegrals source;
    source.electron.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.hole.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    source.combined.assign(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    const auto records = triangleGssAvalancheSourceRecords(
        config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor);
    for (const auto& record : records) {
        const Real electronNodeSource = 0.5 * record.electronSourceIntegral;
        const Real holeNodeSource = 0.5 * record.holeSourceIntegral;
        source.electron[record.node0] += electronNodeSource;
        source.electron[record.node1] += electronNodeSource;
        source.hole[record.node0] += holeNodeSource;
        source.hole[record.node1] += holeNodeSource;
        source.combined[record.node0] += record.node0SourceIntegral;
        source.combined[record.node1] += record.node1SourceIntegral;
    }
    return source;
}

inline std::vector<Real> currentDensityAvalancheSourceIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    return currentDensityAvalancheSourceComponentIntegrals(
        config, impact, mobilityConfig, mobility, edgeCells, mesh, doping,
        cellMaterials, psi, phin, phip, n, p, ni, Vt, fieldFactor).combined;
}

inline std::vector<Real> sgEdgeCurrentAvalancheSourceIntegrals(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& edgeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    const VectorXd&                    psi,
    const VectorXd&                    phin,
    const VectorXd&                    phip,
    const VectorXd&                    n,
    const VectorXd&                    p,
    const std::vector<Real>&           ni,
    Real                               Vt,
    Real                               fieldFactor = 1.0)
{
    std::vector<Real> source(mesh.numNodes(), 0.0);
    const auto records = sgEdgeCurrentAvalancheSourceRecords(
        config,
        impact,
        mobilityConfig,
        mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        psi,
        phin,
        phip,
        n,
        p,
        ni,
        Vt,
        fieldFactor);
    for (const auto& record : records) {
        addMappedEdgeSourceToNodes(
            config, source, edgeCells, mesh, record,
            record.node0SourceIntegral,
            record.node1SourceIntegral,
            record.edgeSourceIntegral);
    }
    return source;
}

inline Real impactIonizationGenerationRate(
    const ImpactIonizationModelConfig& config,
    const ImpactIonizationModel&       impact,
    const MobilityModelConfig&         mobilityConfig,
    const MobilityModel&               mobility,
    const std::vector<std::vector<Index>>& nodeCells,
    const DeviceMesh&                  mesh,
    const DopingModel&                 doping,
    const std::vector<Material>&       cellMaterials,
    Index                              nodeId,
    Real                               electricField,
    Real                               electronDrivingField,
    Real                               holeDrivingField,
    Real                               n,
    Real                               p)
{
    if (config.generation != "current_density")
        return impact.generationRate(electricField, n, p);

    const Real electronImpactField = electronAvalancheDrivingField(
        config, electronDrivingField, electricField, n);
    const Real holeImpactField = holeAvalancheDrivingField(
        config, holeDrivingField, electricField, p);
    const Real alphaN = impact.electronCoefficient(electronImpactField);
    const Real alphaP = impact.holeCoefficient(holeImpactField);
    const Real mun = nodeMobility(
        nodeCells, mesh, doping, mobility, cellMaterials, nodeId, CarrierType::Electron,
        electronImpactField, &mobilityConfig);
    const Real mup = nodeMobility(
        nodeCells, mesh, doping, mobility, cellMaterials, nodeId, CarrierType::Hole,
        holeImpactField, &mobilityConfig);
    return alphaN * mun * std::max(n, 0.0) * std::abs(electronImpactField) +
           alphaP * mup * std::max(p, 0.0) * std::abs(holeImpactField);
}

/// Return average dielectric constant [F/m] for edge @p edgeId.
inline Real edgeEpsilon(const std::vector<std::vector<Index>>& edgeCells,
                        const DeviceMesh&                       mesh,
                        const MaterialDatabase&                 matdb,
                        Index                                   edgeId)
{
    return edgeAvgMaterialProp(edgeCells[edgeId], mesh, matdb,
                               &Material::eps_r, 1.0) * constants::eps0;
}

// ---------------------------------------------------------------------------
// Per-node ni vector
// ---------------------------------------------------------------------------

/// Build per-node intrinsic concentration vector from the material database.
/// For interface nodes shared by multiple regions, uses the first-found value.
inline std::vector<Real> buildNodeNi(const DeviceMesh&       mesh,
                                     const MaterialDatabase& matdb,
                                     Real                    temperature_K = constants::T0)
{
    const Index N = mesh.numNodes();
    std::vector<Real> ni_v(N, 0.0);
    std::vector<bool> found(N, false);
    for (Index c = 0; c < mesh.numCells(); ++c) {
        const auto& cell   = mesh.getCell(c);
        const auto& region = mesh.getRegion(cell.region_id);
        const Real ni_mat = matdb.getMaterial(region.material, temperature_K).ni;
        for (Index nid : cell.node_ids) {
            if (!found[nid]) {
                ni_v[nid]  = ni_mat;
                found[nid] = true;
            }
        }
    }
    return ni_v;
}

/// Validate that the doping model has one entry per mesh node.
inline void validateDopingMeshSize(const DeviceMesh& mesh,
                                   const DopingModel& doping,
                                   const std::string& context)
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            context + ": doping model size does not match mesh node count.");
}

/// Build per-node effective intrinsic concentration including bandgap narrowing.
inline std::vector<Real> buildEffectiveNodeNi(const DeviceMesh&       mesh,
                                              const MaterialDatabase& matdb,
                                              const DopingModel&      doping,
                                              const BandgapNarrowing& bgn,
                                              Real                    thermalVoltage)
{
    const Real temperature_K = thermalVoltage * constants::q / constants::kb;
    std::vector<Real> ni_v = buildNodeNi(mesh, matdb, temperature_K);
    for (Index i = 0; i < mesh.numNodes(); ++i) {
        const Real delta = bgn.deltaEg(doping.totalImpurity(i), 0.0, 0.0);
        ni_v[i] = effectiveIntrinsicDensity(ni_v[i], thermalVoltage, delta);
    }
    return ni_v;
}

/// Validate inputs before building effective intrinsic concentrations.
inline std::vector<Real> buildValidatedEffectiveNodeNi(
    const std::string&             context,
    const DeviceMesh&              mesh,
    const MaterialDatabase&        matdb,
    const DopingModel&             doping,
    const BandgapNarrowingConfig&  bandgapNarrowingConfig,
    Real                           thermalVoltage)
{
    validateDopingMeshSize(mesh, doping, context);
    return buildEffectiveNodeNi(
        mesh,
        matdb,
        doping,
        *makeBandgapNarrowingModel(bandgapNarrowingConfig),
        thermalVoltage);
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

/**
 * @brief Apply strong Dirichlet BCs (row-replacement) to a sparse system.
 *
 * For each constrained node i with prescribed value v:
 *   1. Reduce rhs(k) by A(k,i)*v for all free rows k.
 *   2. Zero column i and row i.
 *   3. Set A(i,i) = 1, rhs(i) = v.
 *
 * The diagonal A(i,i) = 1 is always set explicitly via coeffRef so that
 * the constraint is enforced even if node i had no prior stiffness entries
 * (e.g. insulator nodes in continuity assemblies where all adjacent edges
 * were skipped due to zero mobility).
 */
inline void applyDirichletBC(SparseMatrixd&                         A,
                              VectorXd&                              b,
                              const std::unordered_map<Index, Real>& bcs)
{
    A.makeCompressed();

    // Step 1: propagate prescribed values into free-node RHS
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        for (SparseMatrixd::InnerIterator it(A, i); it; ++it) {
            const int k = static_cast<int>(it.row());
            if (k == i) continue;
            if (bcs.count(static_cast<Index>(k)) == 0)
                b(k) -= it.value() * value;
        }
    }

    // Step 2 & 3: zero Dirichlet rows/cols
    std::unordered_set<int> dirichletSet;
    for (const auto& [nodeId, _] : bcs)
        dirichletSet.insert(static_cast<int>(nodeId));

    for (int col = 0; col < A.outerSize(); ++col) {
        const bool colIsDirichlet = dirichletSet.count(col) > 0;
        for (SparseMatrixd::InnerIterator it(A, col); it; ++it) {
            const int  row          = static_cast<int>(it.row());
            const bool rowIsDirichlet = dirichletSet.count(row) > 0;
            if (rowIsDirichlet || colIsDirichlet)
                it.valueRef() = (row == col) ? 1.0 : 0.0;
        }
    }
    A.prune(0.0);

    // Explicitly ensure A(i,i) = 1 even if the node had no prior entries.
    for (const auto& [nodeId, value] : bcs) {
        const int i = static_cast<int>(nodeId);
        A.coeffRef(i, i) = 1.0;
        b(i) = value;
    }
}

} // namespace vela::detail
