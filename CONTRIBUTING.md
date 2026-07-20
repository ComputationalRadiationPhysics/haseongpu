# Contribution Guidelines

This document is for developers and maintainers. It describes validation before
submitting changes to HASEonGPU; it is not user-facing installation guidance.

Use CI as the CPU baseline, but do not treat it as complete Alpaka 3 backend
coverage. Changes that affect kernels, backend selection, build flags,
packaging, examples, or runtime behavior should also be validated with clean
CUDA and HIP Alpaka backend builds on suitable machines before submission.
Record the backend, compiler, toolkit/runtime version, and non-default CMake
options used for that validation.

Recommended validation flow:

1. Start from a clean working tree.
2. Configure a clean build with `-DHASE_BUILD_RELEASE=ON`,
   `-DHASE_TESTING=ON`, and `-DHASE_ENABLE_PYTHON=ON`.
3. Use `-DHASE_NATIVE_OPTIMIZATIONS=OFF` for redistributable binaries or wheels.
4. Build the tree and run `ctest --output-on-failure`.
5. Install the Python package before running Python tests; `pytest` requires
   the packaged `calcPhiASE` executable and backend-probe library.
6. Validate `import HASEonGPU`, then run `pytest tests`.
7. Run `git diff --check` and pre-commit before merging or tagging.

Pre-commit
----------

This project is set up for use with `pre-commit`. Using it will make your code
conform with most of the easily automatable code style guidelines before it
reaches review.

Install the `pre-commit` executable if it is not already available; see
<https://pre-commit.com> for installation options. Then enable the repository
hooks in your working clone:

```bash
cd /path/to/HASEonGPU-working-clone
python -m pip install pre-commit
pre-commit install --install-hooks
```

The configuration installs both `pre-commit` and `pre-push` hooks. Git will run
the configured checks before each commit and push, and will refuse the action if
one of them fails. Some hooks, such as formatters or whitespace cleanup, may
modify the working tree automatically. In that case, inspect the changes, stage
them, and retry:

```bash
git add -u
git commit
```

To validate the full tree manually, run:

```bash
pre-commit run --all-files
```
