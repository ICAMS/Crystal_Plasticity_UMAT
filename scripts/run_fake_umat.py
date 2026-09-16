#!/usr/bin/env python3
"""Build and run a small local driver that calls the Abaqus UMAT directly."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

tmpdir = Path(tempfile.gettempdir())
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXE = Path(tmpdir / "cp_umat_fake_driver")


LOAD_CASES = {
    "ux": 1,
    "uy": 2,
    "uz": 3,
    "shearxy": 4,
    "shearxz": 5,
    "shearyz": 6,
}

MATERIALS = {
    "aluminum": 1,
    "copper": 2,
    "ferrite": 3,
    "austenite": 4,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile and run a fake Abaqus UMAT material-point driver."
    )
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--points", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--material", choices=sorted(MATERIALS), default="copper")
    parser.add_argument("--load", choices=sorted(LOAD_CASES), default="ux")
    parser.add_argument("--strain", type=float, default=1.0e-3)
    parser.add_argument("--dtime", type=float, default=1.0)
    parser.add_argument("--temp", type=float, default=294.0)
    parser.add_argument("--orientations", choices=["same", "mixed"], default="same")
    parser.add_argument("--order", choices=["natural", "shuffled"], default="natural")
    parser.add_argument("--history", type=Path, help="Write every increment: point, step, pnewdt, stress, column-major tangent, STATEV")
    return parser.parse_args()


def build(args: argparse.Namespace) -> None:
    cmd = [
        args.fc,
        "-fopenmp",
        "-frecursive",
        "-ffixed-line-length-132",
        "-J",
        tmpdir,
        str(ROOT / "src" / "umat.f"),
        str(ROOT / "tests" / "unit" / "umat_fake_driver.f90"),
        "-o",
        str(args.exe),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def run(args: argparse.Namespace) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["OMP_DYNAMIC"] = "FALSE"
    cmd = [
        str(args.exe),
        str(args.steps),
        str(args.points),
        str(MATERIALS[args.material]),
        str(LOAD_CASES[args.load]),
        f"{args.strain:.17g}",
        f"{args.dtime:.17g}",
        f"{args.temp:.17g}",
        str(int(args.orientations == "mixed")),
        str(int(args.order == "shuffled")),
        str(args.history.resolve()) if args.history else "",
    ]
    return subprocess.run(
        cmd,
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )


def main() -> None:
    args = parse_args()
    if not args.no_build:
        build(args)
    result = run(args)
    print(result.stdout, end="")


if __name__ == "__main__":
    main()
