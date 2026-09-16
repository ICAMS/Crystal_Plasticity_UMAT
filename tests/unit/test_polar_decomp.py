"""Analytic polar-factor checks; expected failures document the existing cubic defect."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class PolarDecomposition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory(prefix="cp_polar_") as directory:
            work = Path(directory)
            # Compile the production routine verbatim, without unrelated Abaqus stubs.
            source = (ROOT / "src/other_code.f").read_text()
            start = source.index("      recursive subroutine polar_decomp(")
            end = source.index("\n      END", start) + len("\n      END")
            routine = work / "polar.f"
            routine.write_text(source[start:end] + "\n")
            exe = work / "polar"
            subprocess.run([os.environ.get("FC", "gfortran"),
                            "-ffixed-line-length-132", str(routine),
                            str(ROOT / "tests/unit/polar_standard_cases.f90"),
                            "-o", str(exe)], check=True)
            output = subprocess.run([str(exe)], check=True, text=True,
                                    capture_output=True, timeout=30).stdout
        cls.results = {int(row[0]): row[1:] for row in csv.reader(output.splitlines())}

    def check_case(self, case):
        status, finite, *errors = self.results[case]
        self.assertEqual(status, "0")
        self.assertEqual(finite, "T")
        # Existing single-precision constants limit the repeated-stretch case.
        for label, error in zip(("U", "R", "orthogonality", "reconstruction", "symmetry"), errors):
            self.assertLessEqual(float(error), 1e-8, f"{label} error: {error}")

    def test_identity(self):
        self.check_case(1)

    def test_rigid_rotation(self):
        self.check_case(2)

    def test_uniform_dilation(self):
        self.check_case(3)

    def test_uniaxial_stretch(self):
        self.check_case(4)

    @unittest.expectedFailure
    def test_distinct_stretches_known_cubic_defect(self):
        self.check_case(5)

    @unittest.expectedFailure
    def test_rotated_stretch_known_cubic_defect(self):
        self.check_case(6)

    @unittest.expectedFailure
    def test_simple_shear_known_cubic_defect(self):
        self.check_case(7)
