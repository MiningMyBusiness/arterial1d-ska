# Changelog

## 0.1.0 (2026-11-19)

First public release. DOI: [10.5281/zenodo.24374690](https://doi.org/10.5281/zenodo.24374690).

- One-dimensional blood flow in elastic arterial trees: two-step MacCormack scheme, characteristic
  boundary conditions, junctions with one parent and up to three daughters (total- or
  static-pressure continuity), three-element Windkessel or pure-resistance outlets, prescribed
  inflow.
- Reader for the supported subset of openBF's YAML network format; unsupported features raise
  `UnsupportedFeature`. Bundled ADAN56 network and inflow from openBF (Apache-2.0).
- Wall stiffness from a Young's modulus with a prescribed or empirical wall thickness, or from end
  wave speeds; linear or exponential taper.
- Fixed time step with defaults `cfl = 0.9`, `c_headroom = 1.25`, `dx = 2.5 mm`; the maximum
  Courant number reached is reported (`Result.max_courant`) and a warning is raised above 1.
- Benchmarks against all six cases of Boileau et al. (2015), with the paper's error measures and
  regression tests (the reference data are not distributed; see README).
