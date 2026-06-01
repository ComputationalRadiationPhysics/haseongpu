/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/SerialVersion.hpp>
#include <core/types.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace hase::internal
{
    std::mt19937 rng{5489u}; // or your chosen seed
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    constexpr double SMALL_FOR_LOOPS = 1E-06;

    enum RayLife : bool
    {
        Alive = true,
        Dead = false,
    };

    enum class SurfaceKind
    {
        Face0,
        Face1,
        Face2,
        Up,
        Down
    };

    struct CandidateState
    {
        SurfaceKind kind;
        int forbiddenTag;
    };

    struct HitCandidate
    {
        double length;
        int nextTri;
        int nextCellZ;
        int nextForbidden;
    };

    struct PointM
    {
        double x;
        double y;
        double z;
        std::optional<unsigned> ptIndex = std::nullopt;
        std::optional<unsigned> zIndex = std::nullopt;

        [[nodiscard]] double length() const
        {
            return std::sqrt(x * x + y * y + z * z);
        }

        PointM operator/(double scalar) const
        {
            return PointM{x / scalar, y / scalar, z / scalar};
        }

        PointM operator-(PointM const& other) const
        {
            return PointM{x - other.x, y - other.y, z - other.z};
        }
    };

    struct BaseVersionSerialContext
    {
        std::span<double const> p_in;
        std::span<unsigned const> t_in;
        std::span<double const> beta_v;
        std::span<double const> n_x;
        std::span<double const> n_y;
        std::span<int const> neighbors;
        std::span<float const> surface_new;
        std::span<double const> center_x;
        std::span<double const> center_y;
        std::span<unsigned const> n_p;
        std::span<int const> forbidden;

        RayLife ray_life = Alive;

        unsigned nr_points = 0;
        unsigned N_cells = 0;
        unsigned nr_layers = 0;
        int NumRays = 0;
        double z_mesh = 0.0;
        double sigma_a = 0.0;
        double sigma_e = 0.0;
        double N_tot = 0.0;
    };

    std::array<CandidateState, 5> createStates()
    {
        return {{
            {SurfaceKind::Face0, 0},
            {SurfaceKind::Face1, 1},
            {SurfaceKind::Face2, 2},
            {SurfaceKind::Up, 3},
            {SurfaceKind::Down, 4},
        }};
    }

    int faceOffset(BaseVersionSerialContext const& ctx, SurfaceKind kind)
    {
        switch(kind)
        {
        case SurfaceKind::Face0:
            return 0;
        case SurfaceKind::Face1:
            return static_cast<int>(ctx.N_cells);
        case SurfaceKind::Face2:
            return static_cast<int>(2 * ctx.N_cells);
        default:
            return 0;
        }
    }

    PointM genRandPoint(
        BaseVersionSerialContext const& ctx,
        unsigned t_1,
        unsigned t_2,
        unsigned t_3,
        unsigned layer_i)
    {
        double u = dist(rng);
        double v = dist(rng);

        if((u + v) > 1.0)
        {
            u = 1.0 - u;
            v = 1.0 - v;
        }

        double const w = 1.0 - u - v;

        double const z_rand = (layer_i + dist(rng)) * ctx.z_mesh;

        double const x_rand = ctx.p_in[t_1] * u + ctx.p_in[t_2] * v + ctx.p_in[t_3] * w;

        double const y_rand = ctx.p_in[ctx.nr_points + t_1] * u + ctx.p_in[ctx.nr_points + t_2] * v
                              + ctx.p_in[ctx.nr_points + t_3] * w;

        return PointM{x_rand, y_rand, z_rand};
    }

    std::optional<HitCandidate> tryIntersection(
        BaseVersionSerialContext const& ctx,
        CandidateState const& state,
        PointM const& currentPos,
        PointM const& dir,
        int tri,
        int cell_z,
        double currentBestLength)
    {
        switch(state.kind)
        {
        case SurfaceKind::Face0:
        case SurfaceKind::Face1:
        case SurfaceKind::Face2:
            {
                int const offset = faceOffset(ctx, state.kind);
                int const idx = tri + offset;

                double const nominator
                    = (ctx.n_x[idx] * ctx.p_in[ctx.n_p[idx]] + ctx.n_y[idx] * ctx.p_in[ctx.n_p[idx] + ctx.nr_points])
                      - (ctx.n_x[idx] * currentPos.x + ctx.n_y[idx] * currentPos.y);

                double const denominator = ctx.n_x[idx] * dir.x + ctx.n_y[idx] * dir.y;

                if(denominator == 0.0)
                {
                    return std::nullopt;
                }

                double const length_help = nominator / denominator;

                if(length_help <= 0.0 || length_help > currentBestLength)
                {
                    return std::nullopt;
                }

                return HitCandidate{
                    .length = length_help,
                    .nextTri = ctx.neighbors[idx],
                    .nextCellZ = cell_z,
                    .nextForbidden = ctx.forbidden[idx]};
            }

        case SurfaceKind::Up:
            {
                double const denominator = dir.z;
                if(denominator == 0.0)
                {
                    return std::nullopt;
                }

                double const length_help = ((cell_z + 1) * ctx.z_mesh - currentPos.z) / denominator;

                if(length_help <= 0.0 || length_help > currentBestLength)
                {
                    return std::nullopt;
                }

                return HitCandidate{
                    .length = length_help,
                    .nextTri = tri,
                    .nextCellZ = std::min(cell_z + 1, static_cast<int>(ctx.nr_layers) - 2),
                    .nextForbidden = 4};
            }

        case SurfaceKind::Down:
            {
                double const denominator = dir.z;
                if(denominator == 0.0)
                {
                    return std::nullopt;
                }

                double const length_help = (cell_z * ctx.z_mesh - currentPos.z) / denominator;

                if(length_help <= 0.0 || length_help > currentBestLength)
                {
                    return std::nullopt;
                }

                return HitCandidate{
                    .length = length_help,
                    .nextTri = tri,
                    .nextCellZ = std::max(cell_z - 1, 0),
                    .nextForbidden = 3};
            }
        }

        return std::nullopt;
    }

    double propagation(
        BaseVersionSerialContext& ctx,
        PointM const& rand,
        PointM const& p,
        int t_start,
        int mesh_start,
        int /*N_refl*/)
    {
        PointM currentPos = rand;
        PointM vec = p - rand;

        double const norm = vec.length();
        PointM const dir = vec / norm;

        double distance = norm;
        double const distance_total = norm;
        double gain = 1.0;

        int tri = t_start;
        int cell_z = mesh_start;
        int forb = -1;

        while(true)
        {
            double bestLength = distance;
            std::optional<HitCandidate> bestHit;

            for(auto const& state : createStates())
            {
                if(state.forbiddenTag == forb)
                {
                    continue;
                }

                auto hit = tryIntersection(ctx, state, currentPos, dir, tri, cell_z, bestLength);

                if(hit && hit->length < bestLength)
                {
                    bestLength = hit->length;
                    bestHit = hit;
                }
            }

            double const segmentLength = bestHit ? bestHit->length : distance;

            gain *= std::exp(
                ctx.N_tot * (ctx.beta_v[tri + cell_z * ctx.N_cells] * (ctx.sigma_e + ctx.sigma_a) - ctx.sigma_a)
                * segmentLength);

            distance -= segmentLength;

            currentPos.x += segmentLength * dir.x;
            currentPos.y += segmentLength * dir.y;
            currentPos.z += segmentLength * dir.z;

            if(std::fabs(distance) < SMALL_FOR_LOOPS)
            {
                ctx.ray_life = Dead;
                break;
            }

            if(!bestHit)
            {
                break;
            }

            tri = bestHit->nextTri;
            cell_z = bestHit->nextCellZ;
            forb = bestHit->nextForbidden;
        }

        gain /= distance_total * distance_total;
        return gain;
    }

    void calcImportance(
        BaseVersionSerialContext& ctx,
        PointM const& p,
        std::vector<double>& importance,
        std::vector<int>& N_rays,
        int N_reflections)
    {
        int Rays_dump = 0;
        double sum_phi = 0.0;
        double surf_tot = 0.0;

        for(int i_t = 0; i_t < static_cast<int>(ctx.N_cells); ++i_t)
        {
            for(int i_z = 0; i_z < static_cast<int>(ctx.nr_layers) - 1; ++i_z)
            {
                PointM const startPoint{ctx.center_x[i_t], ctx.center_y[i_t], ctx.z_mesh * (i_z + 0.5)};

                double const prop = propagation(ctx, startPoint, p, i_t, i_z, N_reflections);

                auto const prismIndex = i_t + i_z * ctx.N_cells;

                importance[prismIndex] = ctx.beta_v[prismIndex] * prop;
                sum_phi += importance[prismIndex];
            }

            surf_tot += ctx.surface_new[i_t];
        }

        for(int i_t = 0; i_t < static_cast<int>(ctx.N_cells); ++i_t)
        {
            for(int i_z = 0; i_z < static_cast<int>(ctx.nr_layers) - 1; ++i_z)
            {
                auto const prismIndex = i_t + i_z * ctx.N_cells;

                N_rays[prismIndex] = static_cast<int>(std::floor(importance[prismIndex] / sum_phi * ctx.NumRays));

                Rays_dump += N_rays[prismIndex];
            }
        }

        int const rays_left = ctx.NumRays - Rays_dump;

        if(rays_left > 0)
        {
            for(int i_r = 0; i_r < rays_left; ++i_r)
            {
                auto const rand_t = static_cast<int>(dist(rng) * ctx.N_cells);
                auto const rand_z = static_cast<int>(dist(rng) * (ctx.nr_layers - 1));

                N_rays[rand_t + rand_z * ctx.N_cells]++;
            }
        }

        for(int i_t = 0; i_t < static_cast<int>(ctx.N_cells); ++i_t)
        {
            for(int i_z = 0; i_z < static_cast<int>(ctx.nr_layers) - 1; ++i_z)
            {
                auto const prismIndex = i_t + i_z * ctx.N_cells;

                if(N_rays[prismIndex] > 0)
                {
                    importance[prismIndex]
                        = static_cast<double>(ctx.NumRays) * ctx.surface_new[i_t] * ctx.z_mesh / N_rays[prismIndex];
                }
                else
                {
                    importance[prismIndex] = 0.0;
                }
            }
        }
    }

    std::vector<PointM> createPoints(BaseVersionSerialContext const& ctx)
    {
        std::vector<PointM> points;
        points.reserve(ctx.nr_points * ctx.nr_layers);

        for(unsigned point_i = 0; point_i < ctx.nr_points; ++point_i)
        {
            for(unsigned z = 0; z < ctx.nr_layers; ++z)
            {
                points.emplace_back(
                    PointM{
                        ctx.p_in[point_i],
                        ctx.p_in[ctx.nr_points + point_i],
                        z * ctx.z_mesh,
                        std::make_optional(point_i),
                        std::make_optional(z)});
            }
        }

        return points;
    }

    BaseVersionSerialContext createContext(
        hase::core::ExperimentParameters const& experiment,
        hase::core::HostMesh& mesh)
    {
        auto const hostNumberOfTriangles = mesh.numberOfTriangles;
        auto const raysPerSample = experiment.minRaysPerSample;

        auto const spectralIndex = 0u;
        auto const hostSigmaA = experiment.sigmaA.at(spectralIndex);
        auto const hostSigmaE = experiment.sigmaE.at(spectralIndex);

        return BaseVersionSerialContext{
            .p_in = std::span<double const>{mesh.points.data(), mesh.points.size()},
            .t_in = std::span<unsigned const>{mesh.trianglePointIndices.data(), mesh.trianglePointIndices.size()},
            .beta_v = std::span<double const>{mesh.betaVolume.data(), mesh.betaVolume.size()},

            .n_x = std::span<double const>{mesh.triangleNormalsX.data(), hostNumberOfTriangles},
            .n_y = std::span<double const>{mesh.triangleNormalsY.data(), hostNumberOfTriangles},
            .neighbors = std::span<int const>{mesh.triangleNeighbors.data(), 3 * hostNumberOfTriangles},
            .surface_new = std::span<float const>{mesh.triangleSurfaces.data(), hostNumberOfTriangles},
            .center_x = std::span<double const>{mesh.triangleCenterX.data(), hostNumberOfTriangles},
            .center_y = std::span<double const>{mesh.triangleCenterY.data(), hostNumberOfTriangles},
            .n_p = std::span<unsigned const>{mesh.triangleNormalPoint.data(), 3 * hostNumberOfTriangles},
            .forbidden = std::span<int const>{mesh.forbiddenEdge.data(), 3 * hostNumberOfTriangles},

            .ray_life = Alive,

            .nr_points = mesh.numberOfPoints,
            .N_cells = hostNumberOfTriangles,
            .nr_layers = mesh.numberOfLevels,
            .NumRays = static_cast<int>(raysPerSample),
            .z_mesh = mesh.thickness,
            .sigma_a = hostSigmaA,
            .sigma_e = hostSigmaE,
            .N_tot = mesh.nTot};
    }


} // namespace hase::internal

namespace hase::core
{

    BaseVersionSerial::BaseVersionSerial(ExperimentParameters const& experiment, HostMesh& mesh, Result& result)
        : m_experiment(experiment)
        , m_mesh(mesh)
        , m_result(result)
    {
    }

    void BaseVersionSerial::mainLoop(
        hase::internal::BaseVersionSerialContext& ctx,
        std::vector<hase::internal::PointM> const& points,
        std::vector<double>& dndtAse,
        std::vector<double>& betaCells,
        float hostCrystalFluorescence,
        uint32_t minSampleI,
        uint32_t maxSampleI)
    {
        auto importance = std::vector(ctx.N_cells * (ctx.nr_layers - 1), 0.0);
        auto N_rays = std::vector(ctx.N_cells * (ctx.nr_layers - 1), 0);

        for(hase::internal::PointM const& sample_i : points)
        {
            auto const sampleIndex = sample_i.ptIndex.value() + sample_i.zIndex.value() * ctx.nr_points;

            if(sampleIndex < minSampleI || sampleIndex > maxSampleI)
            {
                continue;
            }
            if(!sample_i.ptIndex.has_value() || !sample_i.zIndex.has_value())
            {
                throw std::runtime_error("Point initialization failed!");
            }

            int realNumRays = 0;
            hase::internal::calcImportance(ctx, sample_i, importance, N_rays, 1);

            for(unsigned triangle_i = 0; triangle_i < ctx.N_cells; ++triangle_i)
            {
                unsigned const t_1 = ctx.t_in[triangle_i];
                unsigned const t_2 = ctx.t_in[ctx.N_cells + triangle_i];
                unsigned const t_3 = ctx.t_in[2 * ctx.N_cells + triangle_i];

                for(unsigned layer_i = 0; layer_i < ctx.nr_layers - 1; ++layer_i)
                {
                    auto const prismIndex = triangle_i + layer_i * ctx.N_cells;

                    realNumRays += N_rays[prismIndex];

                    for(int ray_i = 0; ray_i < N_rays[prismIndex]; ++ray_i)
                    {
                        hase::internal::PointM randPoint = hase::internal::genRandPoint(ctx, t_1, t_2, t_3, layer_i);

                        double const gain = hase::internal::propagation(
                            ctx,
                            randPoint,
                            sample_i,
                            static_cast<int>(triangle_i),
                            static_cast<int>(layer_i),
                            1);


                        // this causes precision problems if phi[sampleIndex] >>> gain * ctx.beta_v[prismIndex] *
                        // importance[prismIndex]
                        m_result.phiAse[sampleIndex] += gain * ctx.beta_v[prismIndex] * importance[prismIndex];
                    }
                }
            }


            if(realNumRays > 0)
            {
                m_result.phiAse[sampleIndex] /= realNumRays;
            }
        }

        for(hase::internal::PointM const& sample_i : points)
        {
            auto const sampleIndex = sample_i.ptIndex.value() + sample_i.zIndex.value() * ctx.nr_points;
            if(sampleIndex < minSampleI || sampleIndex > maxSampleI)
            {
                continue;
            }
            m_result.phiAse[sampleIndex] /= 4.0 * M_PI;

            double const gainPerDensity = betaCells[sampleIndex] * (ctx.sigma_e + ctx.sigma_a) - (ctx.sigma_a);
            m_result.phiAse[sampleIndex] *= (ctx.N_tot / hostCrystalFluorescence);
            m_result.dndtAse[sampleIndex] = gainPerDensity * m_result.phiAse[sampleIndex];
        }

        std::printf("\ncalculations finished, giving back the data\n");
    }

    void BaseVersionSerial::operator()(uint32_t const minSampleI, uint32_t const maxSampleI)
    {
        auto context = hase::internal::createContext(m_experiment, m_mesh);

        auto points = hase::internal::createPoints(context);

        this->mainLoop(
            context,
            points,
            m_result.dndtAse,
            m_mesh.betaCells,
            m_mesh.crystalTFluo,
            minSampleI,
            maxSampleI);
    }

} // namespace hase::core
