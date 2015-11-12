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


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once

// STL
#include <iostream>

// Alpaka
#include <alpaka/alpaka.hpp>

// HASEonGPU
#include <mesh.hpp>
#include <propagate_ray.hpp>
#include <types.hpp>           /* SEED */
#include <RandomGenerator.hpp> /* RandomGenerator */


/**
 * @brief uses a given importance distribution to decide how many rays will be launched from each prism
 *
 * @param *raysDump will contain the number of rays which were mapped to a specific prism
 *
 * for other parameters, see documentation of importanceSampling()
 */
struct DistributeRaysByImportance {
    template <typename T_Acc,
              typename T_Mesh>
    ALPAKA_FN_ACC void operator()( T_Acc const &acc,
                                   T_Mesh const &mesh,
                                   unsigned *raysPerPrism,
                                   double *importance,
                                   float *sumPhi,
                                   const unsigned raysPerSample,
                                   unsigned *raysDump) const {

        auto globalThreadExtent = alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);
        auto globalThreadIdx    = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc);

        auto nThreads             = globalThreadExtent[0];
        unsigned tid              = globalThreadIdx[0];
        unsigned reflection_i     = globalThreadIdx[1];
        unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;
        unsigned raysDumpPart = 0;

        for(unsigned prism_i = tid; prism_i < mesh.numberOfPrisms; prism_i += nThreads){        
            unsigned i = prism_i + reflectionOffset;
            raysPerPrism[i] = (unsigned) floor(importance[i] / (*sumPhi) * raysPerSample);
            raysDumpPart += raysPerPrism[i];
            if(raysPerPrism[prism_i + reflectionOffset] > raysPerSample){
                printf("importance: %f sumPhi: %f raysPerPrism[%d]: %d (max %d)\n",
                       importance[i],
                       *sumPhi, i,
                       raysPerPrism[prism_i+reflectionOffset],
                       raysPerSample);

            }

            assert(raysPerPrism[i] <= raysPerSample);

        }

        alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, raysDump, raysDumpPart);

    }

};

/**
 * @brief takes a number of rays and distributes them randomly over the available prisms
 *        Warning: Does not distribute to reflection slices !!!
 *
 * @param *raysPerPrism the number of rays for each prism (will be changed)
 * @param *raysDump the number of rays which were already distributed
 *
 * for other parameters, see documentation of importanceSampling()
 *
 */
struct DistributeRemainingRaysRandomly {
    template <typename T_Acc,
              typename T_Mesh>
    ALPAKA_FN_ACC void operator()(T_Acc const &acc,
                                  T_Mesh const &mesh,
                                  unsigned *raysPerPrism,
                                  unsigned const raysPerSample,
                                  unsigned const *raysDump) const {

        auto globalThreadExtent = alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);
        auto globalThreadIdx    = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc);
        auto nThreads           = globalThreadExtent.prod();
        auto blockIndex         = alpaka::idx::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0];
        RandomGenerator<T_Acc> rand(acc, SEED, blockIndex);

        assert(raysPerSample >= raysDump[0]);
        unsigned raysLeft = raysPerSample - raysDump[0];

        for(unsigned tid = globalThreadIdx[0]; tid < raysLeft; tid += nThreads){
            int rand_t = static_cast<int>(rand() * mesh.numberOfTriangles);
            int rand_z = static_cast<int>(rand() * (mesh.numberOfLevels-1));
            unsigned randomPrism = rand_t + rand_z * mesh.numberOfTriangles;
            alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, &(raysPerPrism[randomPrism]), 1u);

        }

    }

};

/**
 * @brief corrects the importance to match with the randomly distributed rays
 *
 * @param *raysPerPrism the number of rays to be launced for each prism
 *
 * @param *importance the importance for each prism (will be changed)
 *
 * for other parameters, see documentation of importanceSampling()
 */
struct RecalculateImportance {
    template <typename T_Acc,
              typename T_Mesh>
    ALPAKA_FN_ACC void operator()(T_Acc const &acc,
                                  T_Mesh const &mesh,
                                  unsigned *raysPerPrism,
                                  unsigned raysPerSample,
                                  double *importance) const{

        auto globalThreadExtent = alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);
        auto globalThreadIdx    = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc);
        auto nThreads = globalThreadExtent[0];
        unsigned tid  = globalThreadIdx[0];

        unsigned reflection_i = globalThreadIdx[1];
        unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

        for(unsigned prism_i = tid; prism_i < mesh.numberOfPrisms; prism_i += nThreads){

            int startLevel = prism_i/(mesh.numberOfTriangles);
            int startTriangle = prism_i - (mesh.numberOfTriangles * startLevel);

            if(raysPerPrism[prism_i + reflectionOffset] > 0){
                importance[prism_i + reflectionOffset] = raysPerSample * mesh.getTriangleSurface(startTriangle) / (mesh.surfaceTotal * raysPerPrism[prism_i + reflectionOffset]);
            }
            else{
                importance[prism_i + reflectionOffset] = 0;
            }

        }

    }

};



/**
 * @brief calculates a first estimate on the importance of each prism, based on a single ray started in the center of each prism
 *
 * @param *importance will contain the initial importance for each prism
 * @param *sumPhi will contain the cumulative sum of the importance values
 *
 * For other parameters, see documentation of importanceSampling()
 *
 */
struct PropagateFromTriangleCenter {
    template <typename T_Acc,
              typename T_Mesh>
    ALPAKA_FN_ACC void operator()( T_Acc const &acc,
                                   T_Mesh const &mesh,
                                   double *importance,
                                   const unsigned sample_i,
                                   const double sigmaA,
                                   const double sigmaE) const {

        auto globalThreadExtent = alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);
        auto globalThreadIdx    = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc);

        auto nThreads = globalThreadExtent[0];

        unsigned reflection_i           = globalThreadIdx[1];
        unsigned reflections            = (reflection_i + 1) / 2;
        unsigned reflectionOffset       = reflection_i * mesh.numberOfPrisms;
        ReflectionPlane reflectionPlane = (reflection_i % 2 == 0)? BOTTOM_REFLECTION : TOP_REFLECTION;

        unsigned tid = globalThreadIdx[0];
        for(unsigned prism_i = tid; prism_i < mesh.numberOfPrisms; prism_i += nThreads){

            double gain = 0;

            unsigned startLevel    = prism_i / (mesh.numberOfTriangles);
            unsigned startTriangle = prism_i - (mesh.numberOfTriangles * startLevel);
            Point startPoint       = mesh.getCenterPoint(startTriangle, startLevel);
            Point samplePoint      = mesh.getSamplePoint(sample_i);


            //if(prism_i == 0 && reflections == 6 && reflectionPlane == BOTTOM_REFLECTION)
            gain = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA, sigmaE);

            // if(prism_i == 0)
            //          printf("reflections: %u  prism_i: %u gain: %f\n", reflections, prism_i, gain);


            importance[prism_i + reflectionOffset] = mesh.getBetaVolume(prism_i) * gain;
            if(mesh.getBetaVolume(prism_i) < 0 || gain < 0 || importance[prism_i+reflectionOffset] < 0){
                printf("beta: %f importance: %f gain: %f\n", mesh.getBetaVolume(prism_i), importance[prism_i + reflectionOffset], gain);
            }
        }

    }

};


/**
 * @brief Calculates preImportance which needs only to be done once.
 *
 * @param sample_i         Index of sample point for which importance sampling should be done.
 * @param reflectionSlices 1 + (2 * maxReflections) - Coded information about how many
 *                         reflections a ray should do and on which surface to start.
 * @param dMesh            All information about triangles, points, contants.
 *                         Is located on device memory. See mesh.h for details.
 * @param sigmaA           Absorption value of the ray.
 * @param sigmaE           Emission value of the ray.
 * @param preImportance    Importance based on gain from test rays (will be returned as pointer).
 * @param blockDim         Number of threads per block.
 * @param gridDim          Number of blocks per grid.
 *
 */
template <typename T_Acc, typename T_Workdiv, typename T_Stream, typename T_Mesh>
ALPAKA_FN_HOST void importanceSamplingPropagation( T_Stream &stream,
                                                   T_Workdiv const &workdiv,
                                                   const unsigned sample_i,
                                                   T_Mesh &dMesh,
                                                   const double maxSigmaA,
                                                   const double maxSigmaE,
                                                   double *dPreImportance){

    PropagateFromTriangleCenter propagateFromTriangleCenter;

    auto const exec (alpaka::exec::create<T_Acc>(workdiv, propagateFromTriangleCenter, dMesh, dPreImportance, sample_i, maxSigmaA, maxSigmaE));    
    alpaka::stream::enqueue(stream, exec);

}


/**
 * @brief Calculates importance and ray distribution on prisms.
 *        Based on preImportance and triangle surfaces.
 *
 * @param reflectionSlices   1 + (2 * maxReflections) - Coded information about how many
 *                           reflections a ray should do and on which surface to start.
 * @param dMesh              All information about triangles, points, contants. 
 *                           Is located on device memory. See mesh.h for details.
 * @param raysPerSample      Number of rays that should be distributed on gain medium.
 * @param preImportance      Values calculated from importanceSamplingPropagation.
 * @param importance         Final importance values for this number of rays per sample.
 * @param raysPerPrism       Information about how many rays will start in each prism.
 * @param hSumPhi            Global memory where all threads can sum their phiAse on.
 * @param distributeRandomly In case that not all rays can be distributed by importance
 *                           sampling, the rest will be distributed randomly.
 * @param blockDim           Number of threads per block.
 * @param gridDim            Number of blocks per grid.
 *
 */
template <typename T_Acc, typename T_Workdiv, typename T_Stream, typename T_Mesh>
unsigned importanceSamplingDistribution(T_Stream &stream,
                                        T_Workdiv const &importanceWorkdiv,
                                        T_Workdiv const &workdiv,
                                        const T_Mesh &dMesh,
                                        const unsigned raysPerSample,
                                        double *dPreImportance,
                                        double *dImportance,
                                        unsigned *dRaysPerPrism,
                                        const float sumPhi,
                                        const bool distributeRandomly){

    using DevAcc  = alpaka::dev::Dev<T_Stream>;
    using DevHost = alpaka::dev::DevCpu;

    DevAcc  devAcc  (alpaka::dev::getDev(stream));
    DevHost devHost (alpaka::dev::cpu::getDev());

    auto hSumPhi   ( alpaka::mem::buf::alloc<float,    std::size_t, std::size_t, DevHost>(devHost, static_cast<std::size_t>(1)));
    auto hRaysDump ( alpaka::mem::buf::alloc<unsigned, std::size_t, std::size_t, DevHost>(devHost, static_cast<std::size_t>(1)));
    auto dSumPhi   ( alpaka::mem::buf::alloc<float,    std::size_t, std::size_t, DevAcc> (devAcc, static_cast<std::size_t>(1)));
    auto dRaysDump ( alpaka::mem::buf::alloc<unsigned, std::size_t, std::size_t, DevAcc> (devAcc, static_cast<std::size_t>(1)));

    alpaka::mem::view::getPtrNative(hSumPhi)[0]   = sumPhi;
    alpaka::mem::view::getPtrNative(hRaysDump)[0] = 0;

    alpaka::mem::view::copy(stream, dSumPhi, hSumPhi, static_cast<std::size_t>(1));
    alpaka::mem::view::copy(stream, dRaysDump, hRaysDump, static_cast<std::size_t>(1));

    // Kernel functors
    DistributeRaysByImportance distributeRaysByImportance;
    DistributeRemainingRaysRandomly distributeRemainingRaysRandomly;
    RecalculateImportance recalculateImportance;

    auto const exec0 (alpaka::exec::create<T_Acc>(importanceWorkdiv,
                                                  distributeRaysByImportance,
                                                  dMesh,
                                                  dRaysPerPrism,
                                                  dPreImportance,
                                                  alpaka::mem::view::getPtrNative(dSumPhi),
                                                  raysPerSample,
                                                  alpaka::mem::view::getPtrNative(dRaysDump)));
    alpaka::stream::enqueue(stream, exec0);


    // alpaka::mem::view::copy(stream, hRaysDump, dRaysDump, static_cast<std::size_t>(1));
    // std::cout << alpaka::mem::view::getPtrNative(hRaysDump)[0] << std::endl;


    // Distribute remaining rays randomly if wanted
    if(distributeRandomly){
        auto const exec1 (alpaka::exec::create<T_Acc>(workdiv,
                                                      distributeRemainingRaysRandomly,
                                                      dMesh,
                                                      dRaysPerPrism,
                                                      raysPerSample,
                                                      alpaka::mem::view::getPtrNative(dRaysDump)));
        alpaka::stream::enqueue(stream, exec1);

        // alpaka::mem::view::copy(stream, hRaysDump, dRaysDump, static_cast<std::size_t>(1));
        // std::cout << "RaysLeft:" << raysPerSample - alpaka::mem::view::getPtrNative(hRaysDump)[0] << std::endl;

        alpaka::mem::view::getPtrNative(hRaysDump)[0] = raysPerSample;




    }
    else {
        alpaka::mem::view::copy(stream, hRaysDump, dRaysDump, static_cast<std::size_t>(1));             
    }

    auto const exec2 (alpaka::exec::create<T_Acc>(importanceWorkdiv,
                                                  recalculateImportance,
                                                  dMesh,
                                                  dRaysPerPrism,
                                                  alpaka::mem::view::getPtrNative(hRaysDump)[0],
                                                  dImportance));
    alpaka::stream::enqueue(stream, exec2);

    return alpaka::mem::view::getPtrNative(hRaysDump)[0];
}
