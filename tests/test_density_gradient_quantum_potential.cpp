#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DensityGradientQuantumPotential.h"

#include <cmath>

using namespace vela;

namespace {

DeviceMesh stripMesh()
{
    DeviceMesh mesh;
    mesh.addNode({0, 0.0, 0.0});
    mesh.addNode({1, 1.0e-9, 0.0});
    mesh.addNode({2, 2.0e-9, 0.0});
    mesh.addNode({3, 0.0, 1.0e-9});
    mesh.addNode({4, 1.0e-9, 1.0e-9});
    mesh.addNode({5, 2.0e-9, 1.0e-9});
    mesh.addRegion({0, "silicon", "Silicon", {0, 1, 2, 3}});
    mesh.addCell({0, CellType::Tri3, 0, {0, 1, 4}});
    mesh.addCell({1, CellType::Tri3, 0, {0, 4, 3}});
    mesh.addCell({2, CellType::Tri3, 0, {1, 2, 5}});
    mesh.addCell({3, CellType::Tri3, 0, {1, 5, 4}});
    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("Density-gradient coefficient follows Ancona-Tiersten form",
          "[density_gradient][quantum]")
{
    DensityGradientQuantumPotentialConfig config;
    config.gamma = 3.6;
    config.effectiveMassRatio = 1.0618016171622988;
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    const Real expected = config.gamma * hbar * hbar /
        (6.0 * config.effectiveMassRatio * constants::m0 * constants::q);
    REQUIRE(densityGradientCoefficientVm2(config) ==
            Catch::Approx(expected).epsilon(1.0e-14));
}

TEST_CASE("Uniform density has zero electron quantum potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    const VectorXd density = VectorXd::Constant(6, 1.0e22);
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V.lpNorm<Eigen::Infinity>() < 1.0e-12);
}

TEST_CASE("Carrier confinement raises the electron quantum potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e20, 1.0e22, 1.0e20, 1.0e20, 1.0e22, 1.0e20;
    DensityGradientQuantumPotentialConfig config;
    config.damping = 0.5;
    config.maxIterations = 80;
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V(1) > 0.0);
    REQUIRE(result.potential_V(4) > 0.0);
    REQUIRE(result.potential_V(0) == 0.0);
    REQUIRE(result.potential_V(3) == 0.0);
}

TEST_CASE("Inactive insulator nodes are pinned and decoupled",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density = VectorXd::Constant(6, 1.0e22);
    density(2) = 0.0;
    density(5) = 0.0;
    const std::vector<bool> active = {true, true, false, true, true, false};
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, active, {{0, 0.0}, {3, 0.0}}, constants::Vt_300);
    REQUIRE(result.converged);
    REQUIRE(result.potential_V(2) == 0.0);
    REQUIRE(result.potential_V(5) == 0.0);
}

TEST_CASE("Density-gradient update cap stabilizes a steep imported potential",
          "[density_gradient][quantum]")
{
    const DeviceMesh mesh = stripMesh();
    VectorXd density(6);
    density << 1.0e5, 1.0e27, 1.0e5, 1.0e5, 1.0e27, 1.0e5;
    VectorXd initial(6);
    initial << 0.0, 0.31, -0.31, 0.0, 0.31, -0.31;
    DensityGradientQuantumPotentialConfig config;
    config.maxIterations = 80;
    config.maxUpdate_V = 0.05;
    const auto result = solveElectronDensityGradientPotential(
        mesh, density, std::vector<bool>(6, true), {{0, 0.0}, {3, 0.0}},
        constants::Vt_300, PhysicalUnitSystem::legacySI(), config, initial);
    REQUIRE(result.potential_V.allFinite());
    REQUIRE(std::isfinite(result.residualInfinityNorm));
}
