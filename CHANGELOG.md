# Changelog ICAMS CP-UMAT

All notable changes to this project will be documented here.

## Release  2026R1

### Changed
- Restructured repository layout: `source/` → `src/`, `Web Interface/` → `web/`, `tests/Input/` → `tests/fixtures/`
- Integration tests moved to `tests/integration/`, unit test placeholder added at `tests/unit/`
- Crystal plasticity parameters extracted from [mat-data](https://github.com/ICAMS/mat-data) repository via [mat-data-handler](https://pypi.org/project/mat-data-handler/) in form of include (.inc) files
- All material parameters are imported into Abaqus .inp from .inc file via *Include
- Post-processing scripts moved from `examples/` to `scripts/`
- Added `requirements.txt`, `CONTRIBUTING.md`
- Subroutine for kinematic hardening / backstress evaluation moved to mod_materials.f; subroutine relies only on SDV no dependence on global arrays
- Implemented calculation of equiv. plastic strain based on Uday's PR in GitHub repo
- Made homogenization multi-thread proof
- Calculations of standard material are now also multi-thread proof. All material and simulation parameters are embedded into node-private contexts


### Known Issues (tracked in GitLab)
- Chaboche kinematic hardening: reference parameters lead to unstable Newton-Raphson; using FA/OW values as workaround
- Temperature-dependent parameters not fully implemented; test case `TMF-A3` not yet created
- Superalloy: problems in generated `.inc` file
- Gradient plasticity: throws error in all tested Abaqus versions
- TRIP mechanism: untested

