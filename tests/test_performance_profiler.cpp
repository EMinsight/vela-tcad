#include <catch2/catch_test_macros.hpp>

#include "vela/core/PerformanceProfiler.h"

#include <chrono>
#include <nlohmann/json.hpp>

using namespace vela;

TEST_CASE("PerformanceProfiler aggregates deterministic stages and counters",
          "[performance_profiler]")
{
    PerformanceProfiler profiler({true, "unused.json"});
    profiler.recordStage("newton.jacobian", std::chrono::nanoseconds(10));
    profiler.recordStage("newton.jacobian", std::chrono::nanoseconds(30));
    profiler.increment("newton.updates", 2);
    profiler.observe("linear.rows", 12.0);
    profiler.observe("linear.rows", 18.0);
    profiler.recordNewtonSolve(true, 2, 1.0, 1.0e-9, "abstol",
                               std::chrono::nanoseconds(50));

    const nlohmann::json json = profiler.toJson();
    REQUIRE(json.at("counters").at("newton.updates") == 4);
    REQUIRE(json.at("counters").at("newton.solve_calls") == 1);
    const auto& stage = json.at("stages").at(0);
    REQUIRE(stage.at("name") == "newton.jacobian");
    REQUIRE(stage.at("calls") == 2);
    REQUIRE(stage.at("total_ns") == 40);
    REQUIRE(stage.at("p50_ns") == 10);
    REQUIRE(stage.at("p95_ns") == 30);
    REQUIRE(json.at("observations").at("linear.rows").at("min") == 12.0);
    REQUIRE(json.at("observations").at("linear.rows").at("max") == 18.0);
    REQUIRE(json.at("newton_solves").at(0).at("iterations") == 2);
}

TEST_CASE("Performance profiling config is disabled by default and validates output",
          "[performance_profiler]")
{
    REQUIRE_FALSE(performanceProfilingConfigFromJson(nlohmann::json::object()).enabled);
    const auto config = performanceProfilingConfigFromJson({
        {"performance_profiling", {
            {"enabled", true},
            {"json_file", "profile.json"}
        }}
    });
    REQUIRE(config.enabled);
    REQUIRE(config.jsonFile == "profile.json");
    REQUIRE_THROWS(performanceProfilingConfigFromJson({
        {"performance_profiling", {{"enabled", true}, {"json_file", ""}}}
    }));
}

TEST_CASE("ScopedPerformanceTimer records only inside an active scope",
          "[performance_profiler]")
{
    PerformanceProfiler profiler({true, "unused.json"});
    {
        ScopedPerformanceTimer inactive("inactive");
    }
    {
        ActivePerformanceProfilerScope active(&profiler);
        ScopedPerformanceTimer timer("active");
    }
    const nlohmann::json json = profiler.toJson();
    REQUIRE(json.at("stages").size() == 1);
    REQUIRE(json.at("stages").at(0).at("name") == "active");
}
