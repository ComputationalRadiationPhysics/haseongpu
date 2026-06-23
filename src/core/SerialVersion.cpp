/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/SerialVersion.hpp>

#include <stdexcept>

namespace hase::core
{

    BaseVersionSerial::BaseVersionSerial(ExperimentParameters const& experiment, HostMesh& mesh, Result& result)
        : m_experiment(experiment)
        , m_mesh(mesh)
        , m_result(result)
    {
    }

    void BaseVersionSerial::operator()(uint32_t const minSampleI, uint32_t const maxSampleI)
    {
        (void) m_experiment;
        (void) m_mesh;
        (void) m_result;
        (void) minSampleI;
        (void) maxSampleI;
        throw std::runtime_error("Serial backend is not available for explicit 3D unstructured meshes.");
    }

} // namespace hase::core
