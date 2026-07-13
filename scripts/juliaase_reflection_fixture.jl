#!/usr/bin/env julia

# Generate the single-Tet4 forward/SRM reference used by HASEonGPU's
# JuliaASE regression fixture.  Keep the physical setup in lockstep with
# tests/data/juliaASE/reflection_surface_reference/reference.json.

using Random

const juliaase_root = length(ARGS) == 1 ? ARGS[1] : error("usage: juliaase_reflection_fixture.jl JULIAASE_ROOT")
const ray_count = parse(Int, get(ENV, "HASE_JULIAASE_FIXTURE_RAYS", "4096"))
const max_reflection_passes = parse(Int, get(ENV, "HASE_JULIAASE_FIXTURE_MAX_PASSES", "4"))
include(joinpath(juliaase_root, "src", "ForwardASE.jl"))

const T = ForwardASE.Types
const S = ForwardASE.SRM
const G = ForwardASE.GPUTransfer
const Sim = ForwardASE.Simulation
const Tallies = ForwardASE.Tallies

const points = Float64[
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
]
const connectivity = reshape(Int32[1, 2, 3, 4], 4, 1)
const face_normals = cat(
    reshape(Float32[inv(sqrt(3.0)), inv(sqrt(3.0)), inv(sqrt(3.0))], 3, 1),
    reshape(Float32[-1.0, 0.0, 0.0], 3, 1),
    reshape(Float32[0.0, -1.0, 0.0], 3, 1),
    reshape(Float32[0.0, 0.0, -1.0], 3, 1);
    dims = 2,
)
const mesh = T.TetMesh(
    points,
    Int8[T.TET4],
    connectivity,
    zeros(Int32, 6, 1),
    Int32[1],
    Float32[1.0 / 6.0],
    reshape(Float32[0.25, 0.25, 0.25], 3, 1),
    reshape(face_normals, 3, 4, 1),
    reshape(Float32[sqrt(3.0) / 2.0, 0.5, 0.5, 0.5], 4, 1),
    fill(Int32(-1), 4, 1),
    # HASE's local-face 0 is the z=0 face.  JuliaASE's Tet4 convention numbers
    # that same face 4 (the face opposite vertex 4).
    reshape(Int32[T.BOUND_STOP, T.BOUND_STOP, T.BOUND_STOP, 11], 4, 1),
    fill(NaN, 3, 3, 1),
)

const boundary_faces = T.BoundaryFaceList(
    Int32[1, 1, 1, 1],
    Int8[1, 2, 3, 4],
    Int32[T.BOUND_STOP, T.BOUND_STOP, T.BOUND_STOP, 11],
    face_normals,
    Float32[sqrt(3.0) / 2.0, 0.5, 0.5, 0.5],
    Float32[
        1.0 / 3.0 0.0 1.0 / 3.0 1.0 / 3.0
        1.0 / 3.0 1.0 / 3.0 0.0 1.0 / 3.0
        1.0 / 3.0 1.0 / 3.0 1.0 / 3.0 0.0
    ],
    Int32[0, 0, 0, 1],
)

# HASE's explicit surface reflectivity wins over Fresnel except under total
# internal reflection.  A constant coating has the same behaviour in JuliaASE:
# its reverse-incidence path checks TIR before looking up the coating table.
const coating = T.CoatingTable(
    "hase-surface-reflectivity-0.65",
    Float32[0.0, 90.0],
    Float32[1030.0, 1030.1],
    fill(0.65f0, 2, 2),
    fill(0.65f0, 2, 2),
    fill(0.35f0, 2, 2),
    fill(0.35f0, 2, 2),
)

const beta = 0.18f0
const sigma_absorption = 0.01f0
const sigma_emission = 0.02f0
const gain = beta * (sigma_absorption + sigma_emission) - sigma_absorption
const source_rate_total = Float64(beta) / 6.0
const state = G.init_simulation_state(
    mesh,
    boundary_faces,
    S.init_srm(length(boundary_faces.tet_ind), 64),
    T.TIER3,
    Float32[gain],
    zeros(Float32, 3, 1),
    zeros(Float32, 1),
    fill(beta, 4),
    source_rate_total;
    domain_n = Float32[1.5],
    bfl_outside_n = Float32[1.0, 1.0, 1.0, 1.0],
    coating_tables = T.CoatingTable[coating],
    face_coating_ind = Int32[0, 0, 0, 1],
    track_polarization = false,
)

const run = Sim.run_passes!(
    state,
    ray_count,
    trues(1),
    1.0,
    1030.0f0,
    MersenneTwister(12345);
    epsilon = 1.0e-5,
    max_passes = max_reflection_passes,
    diverge_streak = 3,
    nthreads = 1,
    n_chunks = 1,
)
const phi_ase = Tallies.compute_phi_ase(state, ray_count)

function json_array(values)
    return "[" * join(string.(Float64.(values)), ",") * "]"
end

println("{\"status\":\"", run.status,
        "\",\"passes\":", run.n_passes,
        ",\"rayCount\":", ray_count,
        ",\"phiAse\":", json_array(phi_ase),
        ",\"initialReflectedWeight\":", Float64(sum(run.srm_W_cumulative)),
        ",\"reflectedPassWeightFractions\":", json_array(run.W_fracs),
        "}")
