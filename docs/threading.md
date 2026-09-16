# Threaded UMAT builds

The current UMAT requires full OpenMP compilation and linking, including for
one-CPU reference runs. Selecting `cpus=4` in Abaqus does not enable OpenMP
compilation of the user subroutine. The source still uses `THREADPRIVATE` for
legacy Euler-angle working arrays and the material cache.

`mod_gspt_ini` now validates OpenMP conditional-compilation support. It is called
from both UMAT entry and initial UEXTERNALDB setup. A build that ignores OpenMP
sentinels stops before point initialization with `CP-UMAT requires OpenMP`.
This is an early runtime check, not a compiler/runtime compatibility test.
SIMD-only OpenMP options are not supported, even if they enable sentinels.

## Compiler and linker settings

Keep the Abaqus-provided compiler command, ABI flags, include directories,
shared-library options and libraries. Add the following options to the actual
Fortran compiler command; enable OpenMP on the compiler-driven link command too:

| Compiler | Compile options | Compiler-driven link option |
| --- | --- | --- |
| GNU Fortran (standalone tests) | `-fopenmp -frecursive -ffixed-line-length-132` | `-fopenmp` |
| Intel Fortran, Linux | `-qopenmp -recursive -extend-source 132` | `-qopenmp` |
| Intel Fortran, Windows | `/Qopenmp /recursive /extend-source:132` | `/Qopenmp` |

Use the Intel compiler version supported by your Abaqus installation. Inspect
its effective `compile_fortran` and `link_sl` configuration in the Abaqus
environment and the actual build log. If `link_sl` invokes a linker directly,
configure the compiler's OpenMP runtime libraries there rather than passing a
compiler option to the linker. Do not replace the installation's complete
commands with the table above. Avoid static-local-storage overrides and
OpenMP-disable or SIMD-only options later in the command line. Rebuild the user
library and module files after changing flags; an old library retains its old
thread-storage behavior.

Intel documents that directives become comments without `-qopenmp` or
`/Qopenmp`: [Intel OpenMP support](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-1/add-openmp-support.html).
GNU likewise requires `-fopenmp` for directives and conditional sentinels:
[GNU Fortran OpenMP](https://gcc.gnu.org/onlinedocs/gfortran/OpenMP.html).

## Initialization and ownership

`XI33`, `XI66`, `XI99`, `XInn`, `XI333`, `IB1` and `IB2` are immutable Fortran
parameters with the original values. They need no runtime initialization or
thread-private initialization flag. The `mod_gspt_ini` compatibility entry point
performs no matrix writes and does not reset Euler angles.

The standard solver initializes angles from this call's `PROPS(2:4)` in radians
and loads evolving point matrices from its own `STATEV` into `gp_context`.
Legacy solver paths retain explicit transfers into thread-private module data.
These changes do not modify the polar-decomposition equations or its known
numerical limitations. The previously added F1 assignment remains in place.

## Validation

Run the GNU standalone suite:

```sh
python3 -m unittest discover -s tests/unit -p 'test_*.py' -v
```

The build test requires four actual workers, verifies separate Euler-angle
storage and unchanged tensor constants during repeated initialization, and
checks that compiling without OpenMP produces the intended diagnostic.
The mixed-orientation suite compares serial and threaded stress, tangent and
STATEV histories. The three expected polar-decomposition failures remain known
numerical issues, not successful validations.

Abaqus/Intel validation must still be performed on the production installation:
use the same six-grain mesh and load with identical orientations as a control,
then distinct orientations, on one and four CPUs in thread mode. Compare stress,
state histories, convergence and cutbacks. Repeat for isotropic and kinematic
hardening. GNU tests do not establish that Abaqus workers and the linked Intel
OpenMP runtime provide the intended thread-local storage. This work does not
certify the coupled gradient/TRIP/superalloy paths as thread-safe.
