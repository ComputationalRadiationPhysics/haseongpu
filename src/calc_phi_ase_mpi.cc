/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DISABLE_MPI
    #include <mpi.h>

#include <vector>
#include <iostream>
#include <algorithm>

#include <calc_phi_ase_mpi.hpp>
#include <calc_phi_ase.hpp>
#include <mesh.hpp>
#include <logging.hpp>
#include <progressbar.hpp>
#include <types.hpp>

#include <cuda_runtime_api.h> /* cudaDeviceReset */

// Nodes
#define HEAD_NODE 0

/***
 * Performs ASE computation in MPI mode by statically partitioning the global
 * sample range across all MPI ranks. Each rank (including the HEAD-Node) computes its assigned subset
 * locally, and the partial results are gathered on the head node after all
 * ranks have finished.
 */
float calcPhiAseMPI(
    const ExperimentParameters& experiment,
    const ComputeParameters& compute,
    Mesh& mesh,
    Result& result)
{
    int mpiError = MPI_Init(nullptr, nullptr);
    if (mpiError != MPI_SUCCESS) {
        dout(V_ERROR) << "Error starting MPI program." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, mpiError);
        return 1.0f;
    }

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank != HEAD_NODE) {
        verbosity &= ~V_PROGRESS;
        verbosity &= ~V_STAT;
        verbosity &= ~V_INFO;
    }

    const int totalSamples = static_cast<int>(mesh.numberOfSamples);

    const int base = totalSamples / size;
    const int rem  = totalSamples % size;

    const int localBegin = rank * base + std::min(rank, rem);
    const int localCount = base + (rank < rem ? 1 : 0);
    const int localEnd   = localBegin + localCount;

    float runtime = 0.0f;

    if (localCount > 0) {
        calcPhiAse(
            experiment,
            compute,
            mesh,
            result,
            localBegin,
            localEnd,
            runtime
        );
    }

    std::vector<int> recvCounts;
    std::vector<int> displs;

    if (rank == HEAD_NODE) {
        recvCounts.resize(size);
        displs.resize(size);

        for (int r = 0; r < size; ++r) {
            int rBegin = r * base + std::min(r, rem);
            int rCount = base + (r < rem ? 1 : 0);
            recvCounts[r] = rCount;
            displs[r] = rBegin;
        }
    }

    MPI_Gatherv(
        localCount > 0 ? result.phiAse.data() + localBegin : nullptr,
        localCount,
        MPI_FLOAT,
        rank == HEAD_NODE ? result.phiAse.data() : nullptr,
        rank == HEAD_NODE ? recvCounts.data() : nullptr,
        rank == HEAD_NODE ? displs.data() : nullptr,
        MPI_FLOAT,
        HEAD_NODE,
        MPI_COMM_WORLD
    );

    MPI_Gatherv(
        localCount > 0 ? result.mse.data() + localBegin : nullptr,
        localCount,
        MPI_DOUBLE,
        rank == HEAD_NODE ? result.mse.data() : nullptr,
        rank == HEAD_NODE ? recvCounts.data() : nullptr,
        rank == HEAD_NODE ? displs.data() : nullptr,
        MPI_DOUBLE,
        HEAD_NODE,
        MPI_COMM_WORLD
    );

    MPI_Gatherv(
        localCount > 0 ? result.totalRays.data() + localBegin : nullptr,
        localCount,
        MPI_UNSIGNED,
        rank == HEAD_NODE ? result.totalRays.data() : nullptr,
        rank == HEAD_NODE ? recvCounts.data() : nullptr,
        rank == HEAD_NODE ? displs.data() : nullptr,
        MPI_UNSIGNED,
        HEAD_NODE,
        MPI_COMM_WORLD
    );

    MPI_Gatherv(
    localCount > 0 ? result.dndtAse.data() + localBegin : nullptr,
    localCount,
    MPI_DOUBLE,
    rank == HEAD_NODE ? result.dndtAse.data() : nullptr,
    rank == HEAD_NODE ? recvCounts.data() : nullptr,
    rank == HEAD_NODE ? displs.data() : nullptr,
    MPI_DOUBLE,
    HEAD_NODE,
    MPI_COMM_WORLD
);

    // Gather runtimes to rank 0
    std::vector<float> runtimes;
    if (rank == HEAD_NODE) {
        runtimes.resize(size, 0.0f);
    }

    MPI_Gather(
        &runtime,
        1,
        MPI_FLOAT,
        rank == HEAD_NODE ? runtimes.data() : nullptr,
        1,
        MPI_FLOAT,
        HEAD_NODE,
        MPI_COMM_WORLD
    );

    if (rank == HEAD_NODE) {
        ProgressBar bar;
        for (int i = 0; i < totalSamples; ++i) {
            bar.printFancyProgressBar(mesh.numberOfSamples);
        }
    }
    cudaDeviceReset();
    MPI_Finalize();



    return static_cast<float>(size);
}
#endif
