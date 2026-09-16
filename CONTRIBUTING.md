# Contributing

## Requirements

- Abaqus (v6.14 or later) - needed to run integration tests
- Python 3.8+ with dependencies from `requirements.txt`
- A Fortran compiler (provided by Abaqus; `gfortran` suffices for unit tests)

```bash
pip install -r requirements.txt
```

## Branching

- `main`: stable, tested code only
- Feature branches: `feature-<short-description>`
- Bug fixes: `fix-<short-description>`

Open a merge request against `main`. Link the relevant GitLab issue in the MR description.

## Running the Integration Tests

Each test suite requires the current source to be copied in before running:

```bash
cp -r src tests/integration/test_bkfa/source
cd tests/integration/test_bkfa
abaqus cae noGUI=RunBenchmarkModels.py
```

All test suites must pass before a merge request is approved. 

## Changing the Fortran Source

- `src/umat.f` - entry points only; keep logic out of here
- `src/mod_stress.f` - Newton-Raphson integrator; changes here affect all materials
- `src/mod_wkcoup.f` - weak coupling; changes here affect gradient, TRIP, internal stress, superalloy, and back-stress paths

When modifying a module, re-run the integration tests for every test suite that exercises that code path, not just the one closest to your change.

## Reporting Bugs

Open a GitLab issue with:
- Abaqus version
- Material and flag configuration (`PROPS(5)` value)
- Minimal `.inp` file that reproduces the problem
- Expected vs. actual output
