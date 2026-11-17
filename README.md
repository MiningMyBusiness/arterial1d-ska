# arterial1d

One-dimensional blood flow in trees of elastic arteries, in Python with Numba.

arterial1d is an independent implementation, compatible with the network file format of
[openBF](https://github.com/INSIGNEO/openBF). It solves the same family of equations with a
different numerical scheme (two-step MacCormack finite differences; openBF uses MUSCL finite
volumes). It is not a port of openBF.

> **Status: pre-release.** Validated against one published benchmark so far (ADAN56, below).
> Benchmarks 1–5 of Boileau et al. (2015) are not yet run.

## Install

```bash
pip install -e ".[test]"     # Python >= 3.10; numpy, numba, pyyaml
pytest                       # physics tests; the published-benchmark test skips without its data
```

## Use

```python
import arterial1d as a1

net = a1.load_network(a1.adan56_path(), dx=0.0025)             # any openBF YAML within the supported subset
q_in, period = a1.load_inlet(a1.openbf.inlet_path(a1.adan56_path()))
probes = {"arch": net.node("aortic_arch_I", 0.5), "wrist": net.node("radial_L", 0.95)}
res = net.simulate(q_in, period, n_cycles=15, probes=probes, P_init=net.P_ref)
res.P["wrist"] / a1.MMHG      # pressure [mmHg] at 1 kHz; also res.Q, res.A, per-cycle volumes
```

Networks can also be built directly from `Vessel` objects (radii at both ends, linear or
exponential taper, stiffness from a Young's modulus or from end wave speeds).

## Model

Per vessel, with area A, flow Q, pressure P:

    dA/dt + dQ/dx = 0
    dQ/dt + d(alpha Q^2/A)/dx + (A/rho) dP/dx = -K_R Q/A,     K_R = 2 (gamma + 2) pi mu / rho
    P = P_ext + beta(x) (sqrt(A) - sqrt(A0(x)))

- **Wall:** elastic tube law. From an openBF file, beta follows openBF's convention: the wall
  thickness law of Boileau et al. (2015) and c0 = sqrt(2 E h / (3 rho r)).
- **Friction:** velocity-profile exponent gamma from `gamma_profile` (default 2, i.e. 8 pi mu / rho).
  alpha = 1.
- **Interior scheme:** two-step MacCormack, second order, alternating direction each step. The
  pressure gradient is kept as (A/rho) dP/dx so that a uniform-pressure rest state is preserved
  exactly on tapered vessels.
- **Boundaries:** characteristic compatibility relations in pressure–velocity form, each reduced
  to a scalar Newton solve for the boundary pressure.
  - inlet: prescribed flow Q_in(t);
  - junctions: one parent and up to three daughters; mass conservation and continuity of
    **total** pressure P + rho u²/2 (static-pressure continuity is available, but on ADAN56 it
    gives 4.2% worst-case error against 1.1% for total pressure);
  - outlets: three-element Windkessel R1–C–R2.
- **Time step:** fixed for a run, from the CFL condition at the reference wave speed with headroom
  for its rise with pressure.

## Validation

**ADAN56** (Boileau et al. 2015, benchmark 6; 56 arteries in 77 segments, 31 Windkessel outlets).
Bundled network and inflow from openBF. Mid-vessel pressure over the last of 15 cycles, dx = 2.5 mm,
compared with the published finite-volume (FVM) solution at the 12 published sites. Metric: RMS
difference as % of the reference pulse pressure.

| Solution vs published FVM | median | worst site |
| --- | --- | --- |
| **arterial1d** | **0.51%** | **1.09%** (anterior tibial) |
| published FEM | 0.15% | 0.23% |
| published DCG | 0.12% | 0.18% |
| published LCG | 0.28% | 0.51% |
| published MacCormack | 0.44% | 1.48% |
| published STM | 2.43% | 4.47% |

Diastolic, mean and systolic pressures agree within 0.4 mmHg at every site. Flow RMS differences
are 0.2–1.6% of peak flow. The worst sites are distal: arterial1d's pulse arrives up to 2.7 ms
later than the FVM solution's at the anterior tibial, and this lag does not shrink with grid
(1.25–5 mm) or time-step refinement. Its cause is not yet known.

Reproduce with `python scripts/benchmark_adan56.py --data <folder>` (see *Reference data*).

**Not yet validated:** Boileau et al. benchmarks 1–5 (single pulse, common carotid, upper thoracic
aorta, aortic bifurcation, 37-artery in-vitro network).

## Limitations

- Elastic walls only (no viscoelasticity).
- Trees only: one inlet, one parent per vessel, at most three daughters; no anastomoses or loops.
- Inlet: prescribed flow only (no pressure inlet). Outlets: three-element Windkessel only (no
  reflection-coefficient outlets).
- A subset of openBF's YAML options. The reader raises `UnsupportedFeature` for anything else
  rather than ignoring it.
- Fixed time step per run.
- No speed comparison with other solvers has been made on a matched test. For scale, ADAN56 at
  dx = 5 mm runs in about 0.1 s per cardiac cycle on one laptop core.

## Reference data

The Boileau et al. (2015) supplementary data are the publisher's and are not included here. Download
the supplementary archive from the article page (doi:10.1002/cnm.2732), unpack it, then:

```bash
python scripts/check_reference_data.py /path/to/SuppMaterial_1DBenchmark_Boileau_etal
export ARTERIAL1D_BOILEAU_DIR=/path/to/SuppMaterial_1DBenchmark_Boileau_etal
pytest tests/test_benchmarks.py
```

## Citing

If you use arterial1d, please cite this repository (see `CITATION.cff`) and the sources of the
network and benchmark:

- E. Boileau et al., "A benchmark study of numerical schemes for one-dimensional arterial blood
  flow modelling", Int. J. Numer. Meth. Biomed. Engng 31 (2015) e02732. doi:10.1002/cnm.2732
- P. J. Blanco et al., "An anatomically detailed arterial network model for one-dimensional
  computational hemodynamics", IEEE Trans. Biomed. Eng. 62 (2015) 736–753.
- I. Benemerito, A. Melis, A. Wehenkel, A. Marzo, "openBF: an open-source finite volume 1D blood
  flow solver", Physiol. Meas. 45 (2024) 125002. doi:10.1088/1361-6579/ad9663 (network file format
  and the bundled ADAN56 files)

## Licence

Apache-2.0 (`LICENSE`). The ADAN56 network and inflow files in `src/arterial1d/data/adan56/` are
from openBF, also Apache-2.0; see `NOTICE`.
