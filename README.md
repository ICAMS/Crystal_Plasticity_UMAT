# ICAMS CP-UMAT

Crystal Plasticity User Material (UMAT) subroutine for Abaqus FEM. Implements rate-dependent slip-based plasticity for polycrystalline metals, with support for FCC and BCC crystal structures, multiple isotropic and kinematic hardening models, superalloy mechanisms, strain gradients, and TRIP.


**Organization:** [ICAMS, Ruhr University Bochum](https://www.icams.de), Germany

## Requirements

- [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) v6.14 or later (provides Fortran compiler and Python API)
- Extraction of crystal plasticity parameters for various materials can done with [mat-data-handler](https://pypi.org/project/mat-data-handler/)  
- Current sample materials in the [mat-data](https://github.com/ICAMS/mat-data) repository include: Aluminum, Copper, Ferrite, Austenite, Nickel, Nickel-based Superalloy
- Examples and scripts require Python 3.10+ with [NumPy](https://numpy.org/) and [Matplotlib](https://matplotlib.org/)
- Further optional requirements are [pandas](https://pandas.pydata.org/) and [pylabfea](https://github.com/AHartmaier/pyLabFEA)

All requirements can be installed into your local environment by

```bash
pip install -r requirements.txt
```

## Quick Start

1. Generate a finite element model of the desired geometry as Abaqus .inp file.
2. Prepare the UMAT by copying the contexts of the directory `src` into the working directory and pointing the Abaqus subroutine path for user materials to `umat.f` in your job definition.
3. Extract crystal plasticity parameters as include (.inc) files for used materials from database, e.g. with  

   ```bash
   mat-extract-params copper_generic --outdir includes
   ```  

    **Note:** Upon the first use, you need to initialize the mat-data-handler with  

   ```bash
   mat-extract-params --help
   ```
   
4. The material definition in the .inp file should look like:  

    ```text
    *Material, name=GRAIN1_MAT  
    *Depvar  
        360  
    *User Material, constants=22  
    0., 4.78247643390494, 2.458599547570277, 5.60864095649018, 0., 0., 0., 0.  
    *Include, input="copper_inp_14p.inc"
    ```

for each single crystalline section to which a marterial is assigned, here GRAIN1_MAT. After `*User Material` the number of constant (NPROPS) needs to be specified as the number of parameters in the .inc file for the material plus 8 for the additional elements in the first row. These additional elements contain a legacy index for the specified material (ignored), the three Bunge-Euler angles specifying the crystallographic orientation of the grain in the global coordinate frame in radians, and four spare parameters. Legacy and spare parameter may take any value, but are recommended to be set to "0." After this line of 8 parameters (required by Abaqus), the `*Include` command for the file with the remaining material parameters follows. These parameters will be available as PROPS[9:] within the UMAT.  
**Note:** For downwards compatibility the UMAT accepts a four-parameter legacy input format: a legacy material identifier
and three Bunge-Euler angles in radians. This selects the hard-coded
standard material in `src/mod_material.f`, regardless of the identifier. You must specify the crystal plasticity parameters manually in the Fortran code in this case!  
For either form, the UMAT writes a
material reference to Abaqus unit 6 (the `.dat` file), including the input
constants, all resolved scalar material parameters and flags, and defaults for
inactive blocks. Each distinct PROPS definition is reported once per process,
including when multiple threads initialize it.

5. Submit the job through Abaqus CAE or the command line. Abaqus compiles the Fortran at runtime - no separate build step is needed.  
For command line submission, see `examples/demo_6grain`:

    ```bash
    abaqus job=demo-6gr inp=demo_vox3375_gr6_pbc_strain_geom.inp user=umat.f cpus=4 inter
    ```

For threaded runs and one-CPU references, configure full OpenMP compilation and
linking as described in [Threaded UMAT builds](docs/threading.md). Builds that
ignore OpenMP directives now before point initialization.

## Model generation with Kanapy
To facilitate the setup of microstructure models, it is recommended to use [Kanapy](https://github.com/ICAMS/Kanapy.git), a python package for generating three-dimensional synthetic polycrystals based on characteristic microstructural features developed at ICAMS. The microstructures are built based on statistical information about phase and grain morphologies, given as size distributions and aspect ratio distrubitions of grains and phase regions. Furthermore, crystallographic texture is considered in form of orientation distribution functions (ODF) and misorientation distribution functions (MDF). Kanapy offers tools to analyze EBSD maps with respect to the morphology and texture of microstructures. Based on this experimental data, it generates 3D synthetic microstructures mimicking real ones in a statistical sense.

See [Kanapy/notebooks](https://github.com/ICAMS/Kanapy/tree/master/examples/notebooks) for workflows to automatically generate an Abaqus model referring to the ICAMS CP-UMAP material definitions based on representative volume element (RVE) derived from an EBSD map.

## Running Tests

Integration tests compare simulation output against reference CSV files.

```bash
# Copy current source into the test folder, then run
cp -r src tests/integration/test_bkfa/source
cd tests/integration/test_bkfa
abaqus cae noGUI=RunBenchmarkModels.py
```

Results are written to per-load-case CSV files; `Message_file.txt` contains the pass/fail summary.

The repository also offers unit tests, which are still in the pre-alpha phase.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full test workflow and how to regenerate reference results.

## Repository Layout

```
src/                       Fortran UMAT source
materials/                 Material parameters 
tests/
  integration/             Abaqus regression tests (require Abaqus)
  unit/                    Pure-Fortran unit tests (planned, no Abaqus needed)
  fixtures/                Shared Abaqus input files
examples/                  Example simulation models
scripts/                   Post-processing utilities
docs/                      Technical documentation
```

## Constitutive Model Selection

The options for the constitutive model are selected via an integer selector flag (ISF) in `PROPS(10)`. This flag is ist composed of nibbles (4 bit, max value: 16) in the following table the flag values for each option are given: 

### Isotropic hardening
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 1 |

### Kinematic hardening
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Frederick-Armstrong | 1 | 16 |
| Chaboche | 2 | 32 |
| Ohno-Wang | 3 | 48 |

### Gradient plasticity
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 256 |

### Internal stress for suoeralloy microstructures
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 4 096 |

### Superalloy creep mechanisms
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 65 536 |

### TRIP effect
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 1 048 576 |

### Temperature dependent parameters
| Model | Nibble value | ISF value |
|---|---|---|
| None | 0 | 0 |
| Standard | 1 | 16 777 216 |

The total value of the ISF flag is calculated as the some of the ISF values in the last column of each active option.

ISF = iso\_val + kin\_val + grad\_val + int\_val ...  

See `src/mod_material.f` and `docs/mapping-CP-UMAT-v0_1_0.xlsx` for the full flag layout.


## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Bugs and feature requests go in the [GitHub issue tracker](https://github.com/ICAMS/Crystal_Plasticity_UMAT/issues).

## Authors

Anxin Ma, Martin Boeff, Siwen Gao, Napat Vajragupta, Mahesh Prasad, Alexander Hartmaier

**Institution:** [ICAMS, Ruhr University Bochum](https://www.icams.de), Germany

Contact: <alexander.hartmaier@rub.de>

## License

Source code (including legacy, scripts, examples): [GNU Affero General Public License v3.0](LICENSE)

Documentation (including examples): [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Copyright &copy; 2026 by the Authors
