#pragma once

#include "vela/core/Types.h"
#include "vela/core/UnitScaling.h"
#include "vela/solver/GummelSolver.h"

#include <filesystem>

namespace vela {

DDSolution readDDSolutionStateCsv(
    const std::filesystem::path& path,
    Index expectedNodeCount,
    UnitScalingConfig scaling = UnitScalingConfig{});

void writeDDSolutionStateCsv(
    const std::filesystem::path& path,
    const DDSolution& solution,
    UnitScalingConfig scaling = UnitScalingConfig{});

} // namespace vela
