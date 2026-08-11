# Personal Agent Instructions for HASEonGPU

## Scope and authority

- These instructions apply to this checkout except `juliaASE/` and `thirdParty/`.
- Do not inspect, edit, build, or validate those excluded subtrees unless the user explicitly targets them. When targeted, read and follow their own guidance first.
- Direct user instructions override this file. Do not broaden a request beyond its stated goal.
- Treat `CONTRIBUTING.md` as the validation policy and the relevant files under `docs/source/` as the authority for supported user-facing behavior.

## Start every task from checkout truth

- Read the relevant source, tests, configuration, and documentation before proposing a change.
- Inspect `git status --short --branch` before work. Existing modifications and untracked files belong to the user unless proven otherwise.
- Trace the exact entrypoint named by the user. Do not substitute a similarly named script, wrapper, binary, build tree, or installed package.
- For runtime issues, establish which Python frontend, generated runtime metadata, `calcPhiASE`, C++ openPMD provider, and Python `openpmd_api` are actually loaded. They must come from one compatible runtime.
- Reproduce a reported failure when practical. If reproduction is expensive or unavailable, say which static evidence supports the diagnosis.

## Mandatory explanation and edit gate

- Before editing source, tests, documentation, configuration, or repository instructions, explain the concrete mechanism, the intended edits, and the focused validation.
- Then stop and wait for fresh user approval. This separate approval is required even when the initial request says to implement, fix, or change something.
- A diagnosis, review, explanation, status request, or "do not implement" request never authorizes edits.
- If investigation reveals a materially different fix or expanded scope, explain it and obtain new approval before proceeding.

## Preserve user work

- Never discard, overwrite, reformat, stage, or include unrelated user changes.
- Inspect overlapping diffs before editing a dirty file. If intent cannot be separated safely, stop and ask.
- Apply formatters only to files touched for the approved task. Do not run tree-wide rewriting hooks without explicit approval.
- Avoid destructive Git or filesystem commands. Never use `git reset --hard` or whole-file ours/theirs conflict resolution.
- Resolve only mechanical conflicts autonomously. Surface API, provider, transport, physics, and behavior choices before resolving semantic conflicts.

## Scientific and regression guardrails

- Treat failing tests and numerical differences as diagnostic evidence, not as permission to weaken validation.
- Do not change tests, reference fixtures, expected values, tolerances, RNG seeds, unit conversions, normalization factors, physical constants, or physical models without explaining the scientific consequence and receiving explicit approval for that class of change.
- Preserve deterministic seeds only in the regression paths that require reproducibility; do not leak them into shared runtime defaults.
- Preserve units and array-layout semantics across Python, openPMD, and C++ boundaries.
- For MPI replicated fields, partition global work deterministically, reduce raw integrals, and normalize once after reduction. Do not reduce independently normalized rates.

## Repository conventions

- C++ is C++20 and follows `.clang-format`. Prefer adjacent naming and ownership patterns over introducing new abstractions.
- Python tests use pytest under `tests/python/`; native tests are Catch2 targets generated from `tests/*.cpp` and require `-DHASE_TESTING=ON`.
- Keep repository-facing fixture tooling in Python. Record generator revision, inputs, backend, seed, command/configuration, and checksums when reference data is explicitly regenerated.
- Keep HASE-owned openPMD metadata lower camelCase and update writer, parser, tests, and documentation together when an approved schema change requires it.
- `docs/source/generated/*.rst` is ignored autosummary output. Edit the source `.rst` inventory or prose, never generated pages.
- Keep changes narrowly scoped. Avoid opportunistic cleanup and unrelated renames.

## Build directories

- Never access any directory named `cmake-build-debug`, including for inspection or test discovery.
- Use an existing build only after verifying that its source revision, backend, MPI mode, provider, and relevant CMake options match the task.
- Otherwise use an isolated task-specific directory such as `build/<short-task-name>`.
- A typical test-capable configuration is:

  ```bash
  cmake -S . -B build/<task> -G Ninja \
    -DHASE_BUILD_RELEASE=ON \
    -DHASE_TESTING=ON \
    -DHASE_ENABLE_PYTHON=ON
  ```

- Add only task-required backend, MPI, provider, toolchain, or architecture options. Do not silently change the durable runtime in `build/`.

## Focused validation

- Build and run the smallest target or test selection that exercises the changed behavior. Do not run the full CTest, pytest, pre-commit, or accelerator matrix unless requested.
- Prefer commands of these forms:

  ```bash
  cmake --build build/<task> --target <target>
  ctest --test-dir build/<task> -R '<focused-pattern>' --output-on-failure
  python3 -m pytest <focused-test-path> --tb=short
  pre-commit run --files <touched-files>
  git diff --check
  ```

- Python integration tests require an installed frontend and its packaged runtime artifacts. Verify import and binary/provider resolution before trusting results.
- For documentation changes, validate the edited source with Sphinx warnings treated as errors; do not use committed generated autosummary pages as evidence.
- A CTest result containing `*_NOT_BUILT` is neither a pass nor a regression result. Build the focused target or report it as unverified.
- Report skipped tests, missing hardware, unavailable providers, stale artifacts, and unrun relevant backends explicitly.
- Use CUDA, HIP, MPI, or full-suite validation only when the user requests it or approves it after the focused result shows that broader coverage is needed.

## Rosi and HAL

- For Rosi work, use the `interact-with-rosi` skill. For Alpaka CUDA/HIP/SYCL validation on HAL, use `run-alpaka-backends-on-hal`.
- Use narrow RSYNC snapshots of the local checkout and retrieve only required results. Never use remote Git state or run Git commands on Rosi or HAL.
- Keep remote writes inside the approved workspace, use interactive shells as required by the skill, and record hardware, compiler, runtime/toolkit, backend, revision provenance, and non-default options.

## Git and publication

- Do not stage or commit unless the user explicitly requests a local commit. An explicit commit request authorizes the commit after reviewing its exact diff and scope.
- Ask for immediate confirmation before push, publication, PR/MR creation, issue comments, emails, or any other external mutation.
- Never invoke `gh`. Prepare concise PR, issue, or comment text for the user unless another explicitly approved publication mechanism is requested.
- Preserve requested commit boundaries and ancestry. Do not fold unrelated work into a task commit.

## Handoff

- Lead with the outcome and name the concrete mechanism fixed or diagnosed.
- List the files changed and the focused checks run, including assertion counts or pass/skip/fail results when available.
- State what was not validated and why. Do not imply CPU coverage proves accelerator behavior or that a build-only result proves runtime correctness.
- Mention remaining dirty user changes only to distinguish them from task changes; do not summarize or expose unrelated content.
