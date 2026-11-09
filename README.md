# srtio3-contact-stress

Python calculations supporting Section 5.3.3 of the author's SrTiO3 thesis:
Hertzian spherical indentation and the low-friction scratching approximation
of Lee and Polycarpou (2006).

The script prints contact radius, peak contact pressure, the assumed
subsurface depth, and principal shear-stress estimates. It can also print
worked coefficient calculations. No third-party Python packages are needed.

## Run

Requires Python 3.9 or newer. Open a terminal in the downloaded repository:

```console
python lee_scratch_shear.py --interactive --steps
```

Press Enter to retain a displayed default. On Windows, `py` may be used
instead of `python` if the Python launcher is installed.

Example with explicit parameters:

```console
python lee_scratch_shear.py --load-mn 0.5 --modulus-gpa 225 --radius-nm 2000 --nu 0.24 --mu 0.071 --steps
```

Multiple loads, all evaluated using the same remaining inputs:

```console
python lee_scratch_shear.py --load-mn 0.3 0.5 1.0 --modulus-gpa 225 --radius-nm 2000 --nu 0.24 --mu 0.071
```

Use separate measured onset loads for indentation and scratching. Printing
both estimates at a single load is a model comparison, not a claim that
the two experiments have the same onset load.

## Inputs and Units

| Option | Meaning | Input unit / supported range |
| --- | --- | --- |
| `--load-mn` | Positive normal load(s), P | mN |
| `--modulus-gpa` | Positive reduced contact modulus, E* | GPa |
| `--radius-nm` | Positive effective spherical tip radius, R | nm |
| `--nu` | Specimen Poisson ratio | Dimensionless; 0 <= nu < 0.5 |
| `--mu` | Tangential-to-normal force ratio, Q/P | Dimensionless; 0 <= mu < 0.3 |

E* is the reduced contact modulus, not the specimen Young's modulus.
It already includes elastic compliance of the contacting bodies. The
specimen Poisson ratio is separately needed for the stress distribution.
R is the effective spherical tip radius, not the calculated contact radius a.

Internally, P is converted by multiplying by 1e-3, E* by 1e9, and R by
1e-9 to obtain SI units. Output lengths are nm and stresses are GPa.
The ratios z/a and tau/p0 are dimensionless.

`--material undoped` and `--material nb-doped` select illustrative friction
presets of 0.062 and 0.071. These were estimated from a figure, not computed
from raw friction measurements. They are not intrinsic material constants
or measurements for each tip radius. Replace them with the values and
uncertainty justified for the actual experiments. Explicit `--mu` overrides
the material preset. The default with neither option is mu = 0.071.

## Equations and Numerical Procedure

1. Compute the Hertz contact radius, `a = (3*P*R/(4*E*))^(1/3)`.
2. Compute peak contact pressure, `p0 = 3*P/(2*pi*a^2)`.
3. Use Johnson's Eq. (3.45) for axial and lateral centerline stresses.
   Differentiate their half-difference and solve the stationarity condition
   for `zeta = z/a` by bisection.
4. Evaluate the indentation coefficient as `abs(sigma_z-sigma_r)/(2*p0)`.
5. Evaluate A, B, C and D from Lee and Polycarpou's Eqs. (3)-(4) at this
   same depth, with Q/P identified with mu. Their Eq. (6) gives
   `tau_scratch = p0*D/2`.

The bisection stationarity condition is derived from the cited Hertz
stress expressions; it is not presented as a separately numbered equation
in Johnson. `--steps` shows it and the intermediate numerical values.

### Correspondence with Section 5.3.3

The equations were compared against the supplied Chapter 5 draft dated
2026-10-12. The script reproduces the following expressions and coefficients:

| Thesis equation | Script calculation |
| --- | --- |
| (5.2) | Frictionless indentation: `tau = 0.324047*p0` after rounding, for nu = 0.24 |
| (5.3) | Scratching: `tau = (p0/2)*sqrt((A+B)^2 + 4*mu^2*C^2)` |
| (5.4) | Scratching coefficients 0.324581 and 0.324747 for mu = 0.062 and 0.071, respectively, at nu = 0.24 |

Calculations retain the unrounded coefficients internally. This comparison
checks equations and numerical coefficients, not the underlying experimental
data, the statistical averaging of friction measurements, or plotted values.

## Reproducibility Check

For nu = 0.24:

| Quantity | Value |
| --- | ---: |
| Hertz depth z/a | 0.461009187 |
| Indentation tau/p0 | 0.324047106 |
| Scratching tau/p0, mu = 0.062 | 0.324580925 |
| Scratching tau/p0, mu = 0.071 | 0.324746974 |

For P = 0.5 mN, E* = 225 GPa, R = 2000 nm, nu = 0.24 and mu = 0.071:
a = 149.38 nm, p0 = 10.69857 GPa, z = 68.87 nm,
tau_indentation = 3.466839 GPa and tau_scratching = 3.474327 GPa.
These are example inputs, not a validated experimental dataset.

```console
python lee_scratch_shear.py --reference
python -m unittest discover -s tests -v
```

The reference table also shows nu = 0.30, for which the conventional
indentation coefficient is approximately 0.31. For the SrTiO3 calculations
in this script, nu = 0.24 is the default for both loading modes.

## Scope and Limitations

- Homogeneous, isotropic, linear-elastic half-space and spherical contact.
  Single-crystal elastic anisotropy is not represented explicitly.
- The scratching calculation assumes low friction and retains the
  frictionless Hertz depth, following the approximation discussed by Lee
  and Polycarpou. It is not a full-field sliding-contact maximization.
- The calculated principal shear stress is not the resolved shear stress
  on a particular crystallographic slip system.
- These are onset stress estimates, not simulations of dislocation kinetics,
  crack initiation, fracture toughness, or developed plastic ploughing.
- Changing P, E* or R changes the physical stress and length scales, but
  not z/a in this implementation. Changing nu changes z/a. Friction changes
  the scratching coefficient but not the assumed depth.
- At identical P, E*, R and nu, the implemented friction term can only
  increase or leave unchanged the calculated shear coefficient. It cannot
  by itself explain a lower measured scratching onset stress. Such a
  comparison requires the measured onset loads and the model limitations.
- The code rejects mu >= 0.3; Lee's separate surface-yield expressions are
  not implemented. The cutoff is this model's adopted scope, not a universal
  friction threshold for every material and contact configuration.

## References

1. Johnson, K. L. (1985). *Contact Mechanics*. Cambridge University Press.
   Eq. (3.45), p. 62, supplies the centerline stress expressions.
   https://doi.org/10.1017/CBO9781139171731
2. Lee, K. M., and Polycarpou, A. A. (2006). Shear strength determination
   using the nanoscratch technique and its application to thin solid films.
   *Journal of Materials Research*, 21(9), 2304-2313. Eqs. (1)-(6),
   pp. 2305-2306. https://doi.org/10.1557/JMR.2006.0277

## License

This software is released under the MIT License. See LICENSE for details.

