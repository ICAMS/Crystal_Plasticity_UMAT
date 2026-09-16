Conducting tests for ICAMS CP-UMAT
==================================

Local fake UMAT driver
======================

The script `scripts/run_fake_umat.py` builds and runs a small Fortran driver
that calls `UMAT` directly with Abaqus-like arrays. It is intended for fast
unit/smoke tests without Abaqus, including OpenMP-style concurrent calls:

    python3 scripts/run_fake_umat.py --steps 5 --points 16 --threads 4 --material copper --load ux --strain 1e-4

Mixed-orientation thread regression (requires gfortran with OpenMP):

    python3 -m unittest discover -s tests/unit -p 'test_*.py' -v

The test builds once in an isolated temporary directory and compares all six
stress components, all 36 tangent entries, all 240 STATEV entries, and PNEWDT
at every increment against a serial reference for the same point. It uses
24 points with six distinct Euler-angle triples (radians), plus an identical-
orientation control, under elastic and plastic copper loading. Runs cover
1/2/4 threads, natural/static and shuffled/dynamic scheduling, and two repeats.
Tolerance: max(1e-10, 1e-10 * max(abs(reference), abs(actual))).
Nonfinite values, missing records, and requested cutbacks fail the test.
Compiler errors are failures, not silently skipped tests. Set FC to select a
GNU-compatible Fortran compiler. This test does not emulate Abaqus threading
or validate the production compiler configuration; a passing run cannot rule
out a timing-dependent race.

A deterministic ownership test also interleaves two grain preparations, poisons
legacy module scratch, and checks initial rotations, Fe*Fp, and recovery from
STATEV. The standard solver loads its point context directly; legacy paths
retain a thread-private compatibility boundary. Shared tensor constants are now immutable. See docs/threading.md for build
requirements and the remaining production Abaqus validation.

To inspect a mixed-orientation history manually:

    python3 scripts/run_fake_umat.py --points 24 --threads 4 --orientations mixed --order shuffled --history /tmp/umat-history.csv

History rows have no header: point, increment, PNEWDT, stress(6), tangent(36,
Fortran column-major order), STATEV(240). Normal stdout retains the existing
summary CSV format. The driver stops a point's history at its first cutback;
unexecuted records are zero and rejected by the regression checker.

The driver currently targets standard material behavior. It does not replace
the Abaqus benchmark jobs for coupled gradient, TRIP, superalloy, or output
database validation.

1. Copy files from folder "Input" to case-specific subfolders, i.e. test_int_sup, test_bkow, etc.
2. Create "source" folder with the current UMAT code to be tested against the code in "source-ref" 
3. Make sure all flags are set correctly for the selected test case; compare against settings in "source-ref"
4. Execute "abaqus cae noGUI=RunBenchmarkModels.py" to run test cases. Results will be automatically evaluated an
   copared against results in "*_Ref.csv" files
5. Check Message.txt for maximum errors; check case-specific CSV-files for all errors from this case


Creating reference for new case
===============================

1. Create folder for new test case and copy all files in "Input" to this folder
2. Create a source folder containing the CP-UMAT code for the reference case
3. Make sure all flags are set correctly for the desired test case
4. Edit the file "GenBenchmarks.py" and activate the desired load cases, uncomment others
5. Execute "abaqus cae noGUI=GenBenchmarks.py" to generate the "*_Ref.csv" files containing the reference results
   against which other code versions can be compared
6. Rename "source" into "source-ref"

Polar decomposition standard cases
=================================

    python3 -m unittest discover -s tests/unit -p 'test_polar_decomp.py' -v

The production polar_decomp routine is compiled verbatim and compared against
prescribed analytic factors. Checks cover U, R, R^T R = I, F = R U, symmetry,
finite outputs, and success status. Absolute tolerance is 1e-8 for these
order-one cases because the routine still contains single-precision constants.
Identity, rigid rotation, uniform dilation, and uniaxial stretch pass after
initializing F1 = X1 + sqrt(X2).

Three explicitly expected failures document existing algorithm limitations:
three distinct stretches, rotated nonuniform stretch, and simple shear. They
are NOT validation successes. Clamping a negative cubic discriminant to zero
loses information; fixing F1 alone does not make the general algorithm correct.
Remove the expectedFailure markers when the invariant solver is repaired or
replaced. Unexpected successes will fail the suite to prompt that update.
