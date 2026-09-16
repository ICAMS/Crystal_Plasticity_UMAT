"""Check immutable initialization and refusal of an OpenMP-disabled build."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class ThreadBuild(unittest.TestCase):
    def test_initialization_and_build_guard(self):
        with tempfile.TemporaryDirectory(prefix="cp_thread_build_") as directory:
            work = Path(directory)
            module = work / "module.f"
            module.write_text(
                "      module globalvalue\n"
                "      integer, parameter :: Nslp_mx=60\n"
                "      end module globalvalue\n" +
                (ROOT / "src/mod_gaussp.f").read_text())
            compiler = os.environ.get("FC", "gfortran")
            base = [compiler, "-ffixed-line-length-132", "-J", str(work)]
            exe = work / "enabled"
            subprocess.run(base + ["-fopenmp", str(module),
                           str(ROOT / "tests/unit/thread_initialization.f90"),
                           "-o", str(exe)], check=True, cwd=work)
            subprocess.run([str(exe)], check=True, timeout=30)
            disabled = work / "disabled.f90"
            disabled.write_text("program disabled\nuse mod_gaussp\n"
                                "call mod_gspt_ini()\nend program\n")
            exe = work / "disabled"
            subprocess.run(base + [str(module), str(disabled), "-o", str(exe)],
                           check=True, cwd=work)
            result = subprocess.run([str(exe)], text=True, capture_output=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("CP-UMAT requires OpenMP", result.stdout + result.stderr)
