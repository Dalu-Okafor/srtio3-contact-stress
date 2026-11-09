"""Numerical regression and consistency checks; no third-party packages."""

import math
from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import lee_scratch_shear as model


class ContactStressTests(unittest.TestCase):
    def calculate(self, **changes):
        inputs = dict(load_mn=0.5, reduced_modulus_gpa=225.0,
                      tip_radius_nm=2000.0, poisson_ratio=0.24,
                      friction_coefficient=0.071)
        inputs.update(changes)
        return model.calculate_lee_shear_stress(**inputs)

    def test_reference_coefficients(self):
        r = self.calculate()
        self.assertAlmostEqual(r.depth_ratio, 0.461009186883, places=10)
        self.assertAlmostEqual(r.hertz_tau_over_p0, 0.324047106015, places=10)
        self.assertAlmostEqual(r.lee_tau_over_p0, 0.324746974242, places=10)

    def test_conventional_poisson_ratio(self):
        r = self.calculate(poisson_ratio=0.30)
        self.assertAlmostEqual(r.depth_ratio, 0.480864485, places=8)
        self.assertAlmostEqual(r.hertz_tau_over_p0, 0.310020502, places=8)

    def test_undoped_coefficient(self):
        self.assertAlmostEqual(self.calculate(friction_coefficient=0.062)
                               .lee_tau_over_p0, 0.324580925, places=8)

    def test_zero_friction_matches_hertz(self):
        r = self.calculate(friction_coefficient=0.0)
        self.assertAlmostEqual(r.lee_maximum_shear_gpa,
                               r.hertz_maximum_shear_gpa, places=12)

    def test_hertz_depth_against_independent_grid(self):
        # Direct stress differences on a dense grid, not the derivative solver.
        for nu in (0.0, 0.24, 0.30, 0.49):
            z = model.hertz_maximum_shear_depth_ratio(nu)
            sz, sr = model.hertz_centerline_stresses(nu, z)
            expected = abs(sz - sr) / 2
            grid_max = max(
                abs(x - y) / 2
                for x, y in (model.hertz_centerline_stresses(nu, i / 10000)
                             for i in range(1, 30001))
            )
            self.assertLessEqual(grid_max, expected + 1e-12)
            self.assertLess(expected - grid_max, 1e-8)

    def test_tensor_eigenvalue_difference(self):
        # Independent symmetric 2x2 stress-block eigenvalues at the assumed depth.
        r = self.calculate()
        sz, sr = model.hertz_centerline_stresses(0.24, r.depth_ratio)
        z = r.depth_ratio
        t = 0.071 * (-1 + 1.5*z*math.atan(1/z) - z*z/(2*(1+z*z)))
        midpoint = (sz + sr) / 2
        radius = math.hypot((sr - sz) / 2, t)
        principal = (midpoint - radius, sr, midpoint + radius)
        self.assertAlmostEqual((max(principal)-min(principal))/2,
                               r.lee_tau_over_p0, places=12)

    def test_si_unit_conversion(self):
        r = self.calculate()
        a_m = (3 * 0.0005 * 0.000002 / (4 * 225000000000)) ** (1/3)
        p_pa = 3 * 0.0005 / (2 * math.pi * a_m*a_m)
        self.assertAlmostEqual(r.contact_radius_nm, a_m * 1e9, places=10)
        self.assertAlmostEqual(r.maximum_pressure_gpa, p_pa / 1e9, places=10)

    def test_load_scaling(self):
        a, b = self.calculate(), self.calculate(load_mn=4.0)
        self.assertAlmostEqual(b.contact_radius_nm/a.contact_radius_nm, 2)
        self.assertAlmostEqual(b.maximum_pressure_gpa/a.maximum_pressure_gpa, 2)
        self.assertEqual(a.depth_ratio, b.depth_ratio)

    def test_radius_scaling(self):
        a, b = self.calculate(), self.calculate(tip_radius_nm=16000)
        self.assertAlmostEqual(b.contact_radius_nm/a.contact_radius_nm, 2)
        self.assertAlmostEqual(b.maximum_pressure_gpa/a.maximum_pressure_gpa, 0.25)

    def test_modulus_scaling(self):
        a, b = self.calculate(), self.calculate(reduced_modulus_gpa=1800)
        self.assertAlmostEqual(b.contact_radius_nm/a.contact_radius_nm, 0.5)
        self.assertAlmostEqual(b.maximum_pressure_gpa/a.maximum_pressure_gpa, 4)

    def test_invalid_inputs(self):
        cases = dict(load_mn=(0, -1, float('nan'), float('inf')),
                     reduced_modulus_gpa=(0, -1, float('nan')),
                     tip_radius_nm=(0, -1, float('inf')),
                     poisson_ratio=(-0.1, 0.5, float('nan')),
                     friction_coefficient=(-0.1, 0.3, float('nan')))
        for key, values in cases.items():
            for value in values:
                with self.subTest(parameter=key, value=value):
                    with self.assertRaises(ValueError):
                        self.calculate(**{key: value})

    def test_command_line(self):
        script = str(Path(model.__file__).resolve())
        result = subprocess.run([sys.executable, script, '--load-mn', '0.5', '1',
                                 '--steps', '--reference', '--material', 'undoped'],
                                capture_output=True, text=True, check=True)
        self.assertIn('WORKED COEFFICIENT CALCULATION', result.stdout)
        self.assertIn('REFERENCE COEFFICIENTS', result.stdout)
        self.assertIn('0.324580925', result.stdout)

    def test_explicit_friction_overrides_material(self):
        result = subprocess.run([sys.executable, str(Path(model.__file__).resolve()),
                                 '--material', 'undoped', '--mu', '0'],
                                capture_output=True, text=True, check=True)
        self.assertIn('Friction Q/P = mu  : 0', result.stdout)
        self.assertIn('Friction increase  : 0.000000 %', result.stdout)

    def test_interactive_inputs(self):
        result = subprocess.run([sys.executable, str(Path(model.__file__).resolve()),
                                 '--interactive'], input='0.5\n225\n2000\n0.24\n0.071\n',
                                capture_output=True, text=True, check=True)
        self.assertIn('3.474327', result.stdout)


if __name__ == '__main__':
    unittest.main()
