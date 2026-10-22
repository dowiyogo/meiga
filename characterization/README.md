# Hodoscope Characterization Phase A

This package turns Meiga hodoscope JSON output into geometry-level observables:
angular acceptance, solid angle, geometric efficiency, and the dependence on
the inter-plane distance `D`.

## Modules

- `io.py`: loads Meiga config/output JSON and returns one row per four-plane
  coincidence. Units are mm, rad, ns, and MeV.
- `geometry.py`: stores detector dimensions and maps pixel IDs to centers.
- `acceptance.py`: estimates `T(theta, phi)` and `epsilon(theta)` using the
  EcoMug default `cos^2(theta)` downward muon distribution.
- `solid_angle.py`: compares finite-rectangle analytical solid angle with a
  synthetic Monte Carlo estimate.
- `angular_resolution.py`: provides `sigma_theta ~= d / (sqrt(6) D)` and an
  empirical estimator when truth directions are available.
- `scan_D.py`: scans `D`, writes CSV tables, figures, and `results/REPORT_phaseA.md`.
- `figures_of_merit/run_fom.py`: aggregates existing Phase A/B/C outputs into
  comparison tables and plots without changing upstream physics.

## Normalization

The angular flux convention is

```text
p(theta, phi) dtheta dphi = (3 / 2pi) cos(theta)^2 sin(theta) dtheta dphi
```

for `0 <= theta <= pi/2` and `-pi <= phi < pi`. `pixel_acceptance` returns

```text
T = N_pair_bin / (N_incident_top * P_cos2_bin)
```

so summing `T * P_cos2_bin` over bins and pixel pairs gives the total
four-plane efficiency.

`omega_analytical` returns an area-averaged solid angle in sr. Multiplying it
by the effective top-pixel area gives the two-aperture geometrical factor in
`mm^2 sr`.

## Worked Example

```bash
make characterize
```

or directly:

```bash
python -m characterization.scan_D --d-values 40 80 120 --synthetic-events 5000
```

Outputs are written under `characterization/results/`.

To aggregate the available figures of merit after running the upstream phases:

```bash
make figures-of-merit
```

FoM outputs are written under `data/results/figures_of_merit/`. Missing upstream
artifacts are reported as warnings and remain explicitly unavailable.
