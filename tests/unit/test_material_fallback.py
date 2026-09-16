"""Bounds-checked legacy fallback, explicit materials, and reference output."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class MaterialFallback(unittest.TestCase):
    def test_fallback_and_reporting(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            source = work / 'material.f'
            source.write_text(
                '      module globalvalue\n'
                '      integer, parameter :: Nslp_mx=60\n'
                '      end module\n' +
                (ROOT / 'src/mod_gaussp.f').read_text() +
                (ROOT / 'src/mod_material.f').read_text())
            driver = work / 'driver.f90'
            driver.write_text('''program test
use mod_material
implicit none
type(mat_param_set) :: legacy, explicit
real(8) :: short(4), full(22)
character :: mode
call get_command_argument(1,mode)
short = [99.d0,0.1d0,0.2d0,0.3d0]
full = [2.d0,0.1d0,0.2d0,0.3d0,0.d0,0.d0,0.d0,0.d0, &
        225.d0,1.d0,247000.d0,147000.d0,125000.d0,12.d0, &
        0.001d0,20.d0,20.d0,117.d0,180.d0,1.d0,1.4d0,2.25d0]
if (mode == 's') call extract_material_params(5,full(:5),legacy)
if (mode == 'b') call extract_material_params(17,full(:17),legacy)
! Concurrent first use must produce a single complete reference.
!$omp parallel private(legacy)
call init_material_cached(4,short,legacy)
!$omp end parallel
call init_material_cached(4,short,legacy)
call init_material_cached(22,full,explicit)
if (legacy%ialloy /= 2 .or. legacy%Isf /= 1) stop 2
if (legacy%N_slip /= 12 .or. legacy%spc_grp /= 225) stop 3
if (any(legacy%Mstiff /= explicit%Mstiff)) stop 4
if (any(legacy%HMij /= explicit%HMij)) stop 5
if (any(legacy%IVB_ini /= explicit%IVB_ini)) stop 6
if (legacy%shrt0 /= explicit%shrt0 .or. &
    legacy%pwfl /= explicit%pwfl .or. &
    legacy%pwhd /= explicit%pwhd .or. &
    legacy%hdrt0 /= explicit%hdrt0 .or. &
    legacy%crsss /= explicit%crsss) stop 7
full(11) = 200000.d0
call init_material_cached(22,full,explicit)
if (explicit%c11 /= 200000.d0) stop 8
end program
''')
            exe = work / 'test'
            subprocess.run([os.environ.get('FC', 'gfortran'), '-fopenmp',
                            '-fcheck=all', '-ffixed-line-length-132',
                            str(source), str(driver), '-o', str(exe)],
                           cwd=work, check=True, capture_output=True, text=True)
            result = subprocess.run([str(exe)], check=True, capture_output=True,
                                    text=True, env={**os.environ, 'OMP_NUM_THREADS': '4'})
            self.assertEqual(result.stdout.count('Legacy fallback:'), 1)
            self.assertEqual(result.stdout.count('End UMAT material parameter reference'), 3)
            for name in ['crsss', 'hdrt0', 'pwhd', 'Iwkcoup_temp', 'crss0_slope']:
                self.assertIn(name + ' = ', result.stdout)
            for mode in ['s', 'b']:
                result = subprocess.run([str(exe), mode], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('ERROR:', result.stdout)
