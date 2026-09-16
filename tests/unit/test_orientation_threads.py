"""Run with: python3 -m unittest discover -s tests/unit -p 'test_*.py' -v."""
from __future__ import annotations

import csv
import importlib.util
import math
import os
import subprocess
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("fake_umat", ROOT / "scripts/run_fake_umat.py")
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


class OrientationThreads(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.TemporaryDirectory(prefix="cp_orientation_")
        cls.addClassCleanup(cls.work.cleanup)
        cls.args = SimpleNamespace(
            fc=os.environ.get("FC", "gfortran"),
            exe=Path(cls.work.name) / "driver", steps=5, points=24,
            threads=1, material="copper", load="ux", strain=1e-4,
            dtime=1.0, temp=294.0, orientations="mixed", order="natural",
            history=Path(cls.work.name) / "history.csv",
        )
        # Keep compiler module files isolated from other test/build processes.
        driver.tmpdir = Path(cls.work.name)
        driver.build(cls.args)

    def test_interleaved_orientation_ownership(self):
        work = Path(self.work.name)
        stubs = work / "abaqus_stubs.f90"
        source = (ROOT / "tests/unit/umat_fake_driver.f90").read_text()
        stubs.write_text("subroutine getjobname" + source.split("subroutine getjobname", 1)[1])
        exe = work / "ownership"
        subprocess.run([self.args.fc, "-fopenmp", "-ffixed-line-length-132",
                        "-J", str(work), str(ROOT / "src/umat.f"),
                        str(ROOT / "tests/unit/orientation_ownership.f90"),
                        str(stubs), "-o", str(exe)], cwd=ROOT, check=True)
        subprocess.run([str(exe)], check=True, timeout=30)

    def run_case(self, **changes):
        args = SimpleNamespace(**(vars(self.args) | changes))
        result = driver.run(args)
        records = {}
        with args.history.open() as stream:
            for row in csv.reader(stream):
                key = tuple(map(int, row[:2]))
                values = list(map(float, row[2:]))
                self.assertNotIn(key, records)
                self.assertEqual(len(values), 283)
                self.assertTrue(all(map(math.isfinite, values)), f"nonfinite at {key}")
                self.assertGreaterEqual(values[0], 1.0, f"cutback at {key}\n{result.stdout}")
                records[key] = values
        self.assertEqual(set(records), {(p, s) for p in range(1, args.points+1)
                                       for s in range(1, args.steps+1)})
        return records

    def compare(self, reference, actual):
        for key, expected in reference.items():
            for index, (a, b) in enumerate(zip(expected, actual[key])):
                if index == 0:
                    field = "pnewdt"
                elif index < 7:
                    field = f"stress[{index}]"
                elif index < 43:
                    field = f"tangent[{(index-7)%6+1},{(index-7)//6+1}]"
                else:
                    field = f"STATEV[{index-42}]"
                self.assertTrue(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10),
                                f"point,step={key} {field}: serial={a}, actual={b}")

    def test_orientation_thread_equivalence(self):
        # Fresh process per run also exercises first-use initialization.
        old_schedule = os.environ.get("OMP_SCHEDULE")
        self.addCleanup(self.restore_schedule, old_schedule)
        for strain in (1e-4, 2e-3):  # predominantly elastic and plastic histories
            for orientations in ("same", "mixed"):
                os.environ["OMP_SCHEDULE"] = "static,1"
                base = dict(strain=strain, orientations=orientations)
                reference = self.run_case(**base)
                if orientations == "mixed":
                    stresses = {tuple(reference[p, 5][1:7]) for p in range(1, 7)}
                    self.assertGreater(len(stresses), 1, "orientations must affect response")
                if strain == 2e-3:
                    self.assertTrue(any(v[42+167] > 0 for v in reference.values()),
                                    "plastic case must develop equivalent plastic strain")
                for threads in (1, 2, 4):
                    for order, schedule in (("natural", "static,1"),
                                            ("shuffled", "dynamic,1")):
                        for repeat in range(2):
                            with self.subTest(strain=strain, orientations=orientations,
                                              threads=threads, order=order, repeat=repeat):
                                os.environ["OMP_SCHEDULE"] = schedule
                                self.compare(reference, self.run_case(
                                    **base, threads=threads, order=order))

    @staticmethod
    def restore_schedule(value):
        if value is None:
            os.environ.pop("OMP_SCHEDULE", None)
        else:
            os.environ["OMP_SCHEDULE"] = value


if __name__ == "__main__":
    unittest.main()
