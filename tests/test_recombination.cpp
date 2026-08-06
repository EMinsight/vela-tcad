#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/physics/BandToBandTunnelingModel.h"
#include "vela/equation/AssemblerUtils.h"
#include <nlohmann/json.hpp>

#include <cmath>
#include <stdexcept>
#include <limits>
#include <tuple>
#include <vector>

using namespace vela;

TEST_CASE("Sentaurus E2 band-to-band defaults use the documented SI conversion",
          "[recombination][btbt][sentaurus]")
{
    BandToBandTunnelingConfig config;
    config.model = "e2";
    const BandToBandTunnelingModel model(config);
    const Real field = 1.0e9;
    const Real expected = 3.4e23 * field * field * std::exp(-2.26e9 / field);

    REQUIRE(model.enabled());
    REQUIRE(model.generationRate(field) == Catch::Approx(expected).epsilon(1.0e-13));
    REQUIRE(model.generationRate(-field) == Catch::Approx(expected).epsilon(1.0e-13));
    REQUIRE(model.generationRate(0.0) == 0.0);
}

TEST_CASE("Sentaurus E2 field derivative matches central finite difference",
          "[recombination][btbt][jacobian]")
{
    BandToBandTunnelingConfig config;
    config.model = "e2";
    const BandToBandTunnelingModel model(config);
    const Real field = 8.0e8;
    const Real step = 1.0e3;
    const Real finiteDifference =
        (model.generationRate(field + step) - model.generationRate(field - step)) /
        (2.0 * step);

    REQUIRE(model.generationRateDerivativeField(field) ==
            Catch::Approx(finiteDifference).epsilon(1.0e-9));
}

TEST_CASE("Band-to-band JSON accepts SI and Sentaurus-native parameter units",
          "[recombination][btbt][json]")
{
    const nlohmann::json json = {
        {"model", "e2"},
        {"A_cm_inv_s_inv_V_inv2", 3.4e21},
        {"B_V_per_cm", 22.6e6},
        {"source_integration", "transport_node_lumped"},
        {"jacobian", "frozen_field"},
    };
    const BandToBandTunnelingConfig config =
        bandToBandTunnelingConfigFromJson(json, "test");

    REQUIRE(config.prefactorA_SI == Catch::Approx(3.4e23));
    REQUIRE(config.exponentialB_V_per_m == Catch::Approx(2.26e9));
    REQUIRE(config.sourceIntegration == "transport_node_lumped");
    REQUIRE(config.jacobian == "frozen_field");
    REQUIRE_NOTHROW(BandToBandTunnelingModel(config));
}

TEST_CASE("E2 source integration excludes insulator cells at shared nodes",
          "[recombination][btbt][mixed-material]")
{
    DeviceMesh mesh;
    for (const auto [id, x, y] : std::vector<std::tuple<Index, Real, Real>>{
             {0, 0.0, 0.0}, {1, 1.0, 0.0},
             {2, 0.0, 1.0}, {3, 1.0, 1.0}}) {
        Node node;
        node.id = id;
        node.x = x;
        node.y = y;
        mesh.addNode(node);
    }
    Cell silicon;
    silicon.id = 0;
    silicon.type = CellType::Tri3;
    silicon.region_id = 0;
    silicon.node_ids = {0, 1, 2};
    mesh.addCell(silicon);
    Cell oxide;
    oxide.id = 1;
    oxide.type = CellType::Tri3;
    oxide.region_id = 1;
    oxide.node_ids = {1, 3, 2};
    mesh.addCell(oxide);

    Material si;
    si.name = "Si";
    si.ni = 1.0;
    si.mun = 1.0;
    Material sio2;
    sio2.name = "SiO2";
    const std::vector<Material> cellMaterials{si, sio2};

    VectorXd psi(4);
    psi << 0.0, 1.0, 0.0, 100.0;
    BandToBandTunnelingConfig config;
    config.model = "e2";
    config.prefactorA_SI = 1.0;
    config.exponentialB_V_per_m = 1.0;
    const BandToBandTunnelingModel model(config);
    const std::vector<Real> sources =
        detail::bandToBandGenerationNodeSourceIntegrals(
            model, PhysicalUnitSystem::legacySI(), mesh, cellMaterials, psi, 1.0);
    const Real expectedLumped = std::exp(-1.0) * 0.5 / 3.0;

    REQUIRE(sources[0] == Catch::Approx(expectedLumped));
    REQUIRE(sources[1] == Catch::Approx(expectedLumped));
    REQUIRE(sources[2] == Catch::Approx(expectedLumped));
    REQUIRE(sources[3] == 0.0);

    config.sourceIntegration = "transport_node_lumped";
    const BandToBandTunnelingModel nodalModel(config);
    const std::vector<Real> nodalSources =
        detail::bandToBandGenerationNodeSourceIntegrals(
            nodalModel, PhysicalUnitSystem::legacySI(), mesh,
            cellMaterials, psi, 1.0);
    REQUIRE(nodalSources[0] == Catch::Approx(expectedLumped));
    REQUIRE(nodalSources[1] == Catch::Approx(expectedLumped));
    REQUIRE(nodalSources[2] == Catch::Approx(expectedLumped));
    REQUIRE(nodalSources[3] == 0.0);
}

TEST_CASE("Band-to-band source integration rejects unknown spatial recovery",
          "[recombination][btbt][json]")
{
    BandToBandTunnelingConfig config;
    config.model = "e2";
    config.sourceIntegration = "unknown";
    REQUIRE_THROWS_AS(BandToBandTunnelingModel(config), std::invalid_argument);
}

TEST_CASE("SRH recombination is near zero when n*p equals ni squared", "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"srh"}, 1.0e-7, 2.0e-7));
    const Real ni = 1.0e16;
    const Real n = 2.0e21;
    const Real p = ni * ni / n;

    REQUIRE(model.srhRate(n, p, ni) == Catch::Approx(0.0).margin(1.0e6));
    REQUIRE(model.totalRate(n, p, ni) == Catch::Approx(0.0).margin(1.0e6));
}

TEST_CASE("Auger recombination increases at high carrier concentration", "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"auger"}));
    const Real ni = 1.0e16;

    const Real low = model.augerRate(1.0e21, 1.0e21, ni);
    const Real high = model.augerRate(1.0e24, 1.0e24, ni);

    REQUIRE(low > 0.0);
    REQUIRE(high > low);
}

TEST_CASE("Auger recombination rejects negative coefficients", "[recombination]")
{
    RecombinationModelConfig cfg = recombinationModelConfig({"auger"});
    cfg.augerCn = -1.0e-43;
    REQUIRE_THROWS_AS(RecombinationModel(cfg), std::invalid_argument);

    cfg = recombinationModelConfig({"auger"});
    cfg.augerCp = -1.0e-43;
    REQUIRE_THROWS_AS(RecombinationModel(cfg), std::invalid_argument);
}

TEST_CASE("Auger linearization remains finite for extreme initializer carriers",
          "[recombination]")
{
    RecombinationModel model(recombinationModelConfig({"auger"}));
    const Real ni = 1.0e16;
    const Real n = 1.0e234;
    const Real p = 1.0e234;

    const RecombinationLinearization electron = model.electronLinearization(n, p, ni);
    const RecombinationLinearization hole = model.holeLinearization(n, p, ni);

    REQUIRE(std::isfinite(electron.diagonal));
    REQUIRE(std::isfinite(electron.rhs));
    REQUIRE(std::isfinite(hole.diagonal));
    REQUIRE(std::isfinite(hole.rhs));
}

TEST_CASE("Total recombination is SRH plus Auger", "[recombination]")
{
    RecombinationModel total(recombinationModelConfig({"srh", "auger"}));
    const Real n = 1.0e22;
    const Real p = 2.0e22;
    const Real ni = 1.0e16;

    REQUIRE(total.totalRate(n, p, ni) ==
            Catch::Approx(total.srhRate(n, p, ni) + total.augerRate(n, p, ni)));
}

TEST_CASE("Default recombination parameters match Sentaurus 2018 silicon at 300 K",
          "[recombination][sentaurus]")
{
    const RecombinationModelConfig cfg = recombinationModelConfig({"srh", "auger"});

    REQUIRE(cfg.taun == Catch::Approx(1.0e-5));
    REQUIRE(cfg.taup == Catch::Approx(3.0e-6));
    REQUIRE(cfg.augerCn == Catch::Approx(2.90e-43).epsilon(1.0e-12));
    REQUIRE(cfg.augerCp == Catch::Approx(1.028e-43).epsilon(1.0e-12));
}

TEST_CASE("Default bandgap narrowing interface returns zero", "[bgn]")
{
    NoBandgapNarrowing bgn;
    REQUIRE(bgn.deltaEg(1.0e25, 1.0e24, 1.0e20) == Catch::Approx(0.0));
}

TEST_CASE("Slotboom bandgap narrowing grows effective intrinsic density", "[bgn]")
{
    BandgapNarrowingConfig cfg;
    cfg.model = "slotboom";
    SlotboomBandgapNarrowing bgn(cfg);

    const Real low = bgn.deltaEg(1.0e21, 0.0, 0.0);
    const Real high = bgn.deltaEg(1.0e25, 0.0, 0.0);

    REQUIRE(low >= 0.0);
    REQUIRE(high > low);
    REQUIRE(high == Catch::Approx(0.0833788).epsilon(1.0e-5));

    const Real ni = 1.0e16;
    const Real Vt = 0.025852;
    const Real niEff = effectiveIntrinsicDensity(ni, Vt, high);
    REQUIRE(niEff > ni);
    REQUIRE(niEff == Catch::Approx(ni * std::exp(high / (2.0 * Vt))));
}

TEST_CASE("Effective intrinsic density caps overflow", "[bgn]")
{
    const Real ni = 1.0e16;
    const Real Vt = 0.025852;
    const Real niEff = effectiveIntrinsicDensity(ni, Vt, 1.0e6);

    REQUIRE(std::isfinite(niEff));
    REQUIRE(niEff == std::numeric_limits<Real>::max());
}

TEST_CASE("Bandgap narrowing factory validates model names", "[bgn]")
{
    REQUIRE(makeBandgapNarrowingModel(bandgapNarrowingConfig("none"))->deltaEg(1.0e25, 0.0, 0.0)
            == Catch::Approx(0.0));
    REQUIRE(makeBandgapNarrowingModel(bandgapNarrowingConfig("slotboom"))->deltaEg(1.0e25, 0.0, 0.0)
            > 0.0);
    REQUIRE(makeBandgapNarrowingModel(bandgapNarrowingConfig("old_slotboom"))->deltaEg(1.0e24, 0.0, 0.0)
            == Catch::Approx(0.0424016824523).epsilon(1.0e-10));
    REQUIRE_THROWS_AS(makeBandgapNarrowingModel(bandgapNarrowingConfig("unknown")),
                      std::invalid_argument);
}

TEST_CASE("OldSlotboom matches Sentaurus split ni and positive BGN semantics", "[bgn][sentaurus]")
{
    const Real materialNiFromBandgap = 1.4638914958767616e16;
    const Real dopingAtNref = 1.0e23;
    const Real Vt = 0.025851999786435;

    const auto bgn = makeBandgapNarrowingModel(bandgapNarrowingConfig("old_slotboom"));
    const Real deltaEg = bgn->deltaEg(dopingAtNref, 0.0, 0.0);
    const Real niEff = effectiveIntrinsicDensity(materialNiFromBandgap, Vt, deltaEg);

    REQUIRE(deltaEg == Catch::Approx(0.006363961030678928).epsilon(1.0e-12));
    REQUIRE(niEff / 1.0e6 == Catch::Approx(1.6556319846864e10).epsilon(1.0e-10));
}

TEST_CASE("OldSlotboom reference doping follows unit scaling concentration units", "[bgn][scaling]")
{
    const UnitScalingConfig scaling{UnitScalingMode::UnitScaling};
    const BandgapNarrowingConfig config =
        bandgapNarrowingConfig("old_slotboom", scaling);

    REQUIRE(config.referenceDoping == Catch::Approx(1.0e17));

    const auto bgn = makeBandgapNarrowingModel(config);
    const Real deltaEg = bgn->deltaEg(1.0e17, 0.0, 0.0);
    REQUIRE(deltaEg == Catch::Approx(0.006363961030678928).epsilon(1.0e-12));

    const Real materialNi = 1.4638914958767616e10;
    const Real niEff = effectiveIntrinsicDensity(
        materialNi, 0.025851999786435, deltaEg);
    REQUIRE(niEff == Catch::Approx(1.6556319846864e10).epsilon(1.0e-10));
    REQUIRE(niEff * niEff / 1.0e17 ==
            Catch::Approx(2741.1172687166).epsilon(1.0e-10));
}

TEST_CASE("Fermi statistics adds the Sentaurus apparent BGN correction", "[bgn][fermi][sentaurus]")
{
    const Real correction = fermiStatisticsBandgapCorrection(
        5.1e26, 1.0e24, 2.8e25, 1.04e25, 0.025851999786435);

    REQUIRE(correction == Catch::Approx(0.13920003585558116).epsilon(2.0e-12));
    REQUIRE(fermiStatisticsBandgapCorrection(
        0.0, 0.0, 2.8e25, 1.04e25, 0.025851999786435) == Catch::Approx(0.0));
}

TEST_CASE("CarrierStatistics intrinsic density uses temperature_K material path", "[temperature]")
{
    MaterialDatabase matdb;
    const Material& si = matdb.getMaterial("Si");
    const Real ni300 = intrinsicDensity(si, 300.0);
    const Real ni450 = intrinsicDensity(si, 450.0);
    REQUIRE(ni300 == Catch::Approx(si.ni));
    REQUIRE(ni450 > ni300);
}
