"""Comparison with the published solutions of Boileau et al. 2015 (doi:10.1002/cnm.2732).

The reference data are the paper's supplementary material and are not distributed with this
package. Download them yourself (see README), unpack, and point ARTERIAL1D_BOILEAU_DIR (or the
`root` argument) at the folder that contains Benchmark1_SinglePulseModel ... Benchmark6_ADAN56.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from . import openbf
from .fourier import read_inflow_formula
from .metrics import all_errors, on_times
from .solver import MMHG, Network, Vessel

ENV = "ARTERIAL1D_BOILEAU_DIR"

# Column order of Benchmark6_ADAN56/NumData/*_P.txt and *_Q.txt, as ADAN56 labels in openBF's file
ADAN56_SITES = ["aortic_arch_I", "thoracic_aorta_III", "abdominal_aorta_V", "common_carotid_R", "renal_R",
                "common_iliac_R", "internal_carotid_R", "radial_R", "internal_iliac_R",
                "posterior_interosseous_R", "femoral_R_II", "anterior_tibial_R"]
SCHEMES = ["FVM", "FEM", "DCG", "LCG", "McC", "STM"]


def reference_root(root: Optional[os.PathLike] = None) -> Optional[Path]:
    """The unpacked supplementary folder, or None if it is not available."""
    r = Path(root) if root is not None else (Path(os.environ[ENV]) if ENV in os.environ else None)
    if r is None or not (r / "Benchmark6_ADAN56" / "NumData" / "FVM_P.txt").exists():
        return None
    return r


def _last_cycle(path: Path, period: float, scale: float):
    """The last `period` of a published solution, as (time in cycle, values). Time in cycle is
    t mod period: the published runs start at t = 0 at the start of a cycle."""
    d = np.loadtxt(path)
    dt = d[1, 0] - d[0, 0]
    n = int(round(period / dt))
    d = d[-n:]
    t = np.mod(d[:, 0] + 1e-9, period) - 1e-9
    o = np.argsort(t)
    return t[o], d[o, 1:] * scale


def adan56(root: Optional[os.PathLike] = None, *, dx: float = 0.0025, cycles: int = 15,
           junction: str = "total", reference: str = "FVM") -> Dict[str, dict]:
    """Run ADAN56 with the benchmark inflow and compare mid-vessel pressure and flow over the last
    cycle with a published scheme (FVM by default, the reference used in the paper).

    Returns, per site: pressure RMS difference as a fraction of the reference pulse pressure,
    flow RMS difference as a fraction of the reference peak |flow|, and (min, mean, max)
    pressure in mmHg for both solutions.
    """
    r = reference_root(root)
    if r is None:
        raise FileNotFoundError(f"Boileau et al. 2015 supplementary data not found; set {ENV}")
    cfg = openbf.inlet_path(_adan56_yaml())
    net = openbf.load_network(_adan56_yaml(), dx=dx, junction=junction)
    q, period = openbf.load_inlet(cfg)
    res = net.simulate(q, period, cycles, {a: net.node(a, 0.5) for a in ADAN56_SITES}, P_init=net.P_ref)
    if not res.finite:
        raise RuntimeError("simulation did not stay finite")
    n = int(round(period * 1000))                     # simulate() records at 1 kHz, at step ends
    t = np.mod(res.t[-n:] + 1e-9, period) - 1e-9
    o = np.argsort(t)
    nd = r / "Benchmark6_ADAN56" / "NumData"
    tp, P = _last_cycle(nd / f"{reference}_P.txt", period, 1e3)        # kPa -> Pa
    tq, Q = _last_cycle(nd / f"{reference}_Q.txt", period, 1e-6)       # mL/s -> m^3/s
    tg = np.arange(0.0, period, 0.002)
    out = {}
    for i, a in enumerate(ADAN56_SITES):
        p = np.interp(tg, t[o], res.P[a][-n:][o], period=period)
        f = np.interp(tg, t[o], res.Q[a][-n:][o], period=period)
        pr = np.interp(tg, tp, P[:, i], period=period)
        qr = np.interp(tg, tq, Q[:, i], period=period)
        out[a] = {"P_rms_over_PP": float(np.sqrt(np.mean((p - pr) ** 2)) / (pr.max() - pr.min())),
                  "Q_rms_over_peak": float(np.sqrt(np.mean((f - qr) ** 2)) / np.abs(qr).max()),
                  "P_mmHg": tuple(float(x) / MMHG for x in (p.min(), p.mean(), p.max())),
                  "P_ref_mmHg": tuple(float(x) / MMHG for x in (pr.min(), pr.mean(), pr.max()))}
    return out


def _adan56_yaml():
    from . import adan56_path
    return adan56_path()


# ---- case 2: common carotid artery (Boileau et al. 2015, section 4.2, Table III) -----------------
# Elastic tube, P = Pd + (4/3) sqrt(pi) E h / Ad (sqrt(A) - sqrt(Ad)); three-element Windkessel.
CAROTID = dict(L=0.126, rd=3.0e-3, h=0.3e-3, E=700.0e3, Pd=10.933e3, rho=1060.0, mu=4.0e-3, zeta=2,
               R1=2.4875e8, R2=1.8697e9, C=1.7529e-10)
CAROTID_LAST_CYCLE = 8.8          # README: last cycle 8.8-9.9 s in the 1-D files
CAROTID_SHIFT_3D = 0.056          # the 3-D cycle starts 0.055 s later in the inflow (README), plus 1 ms
                                  # time-label offset; reproduces Table IV exactly for FVM and FEM


def _single_vessel(p, dx, junction=True):
    v = Vessel("vessel", p["L"], p["rd"], p["rd"], None, E=p["E"], h=p["h"], R1=p["R1"], R2=p["R2"], C=p["C"])
    return Network([v], rho=p["rho"], mu=p["mu"], alpha=1.0, friction_coeff_pi=2.0 * (p["zeta"] + 2), dx=dx,
                   P_ven=0.0, P_ref=p["Pd"])


def _published_cycle(path, t0, period):
    """One cycle of a published 1-D solution: (t in cycle, P [Pa], Q [m^3/s], dr [m], dP [Pa])."""
    d = np.loadtxt(path)
    d = d[(d[:, 0] >= t0 - 1e-9) & (d[:, 0] < t0 + period - 1e-9)]
    return d[:, 0] - t0, d[:, 1] * 1e3, d[:, 2] * 1e-6, d[:, 3] * 1e-3, d[:, 4] * 1e3


def carotid(root: Optional[os.PathLike] = None, *, dx: float = 0.0025, cycles: int = 10) -> dict:
    """Run case 2 and compare mid-vessel P, Q, dr and the inlet-outlet pressure drop with the 3-D
    solution and with the six published 1-D schemes, using the paper's error measures (%)."""
    r = reference_root(root)
    if r is None:
        raise FileNotFoundError(f"Boileau et al. 2015 supplementary data not found; set {ENV}")
    base = r / "Benchmark2_CommonCarotidArtery"
    f, period, _ = read_inflow_formula(base / "inflow.txt")
    p = CAROTID
    net = _single_vessel(p, dx)
    probes = {"mid": net.node("vessel", 0.5), "in": net.node("vessel", 0.0), "out": net.node("vessel", 1.0)}
    res = net.simulate(lambda t: f(t) * 1e-6, period, cycles, probes, P_init=p["Pd"])
    if not res.finite:
        raise RuntimeError("simulation did not stay finite")
    n = int(round(period * 1000))
    t = res.t[-n:] - (cycles - 1) * period                  # time in the last cycle (true time)
    mine = {"P": res.P["mid"][-n:], "Q": res.Q["mid"][-n:],
            "dr": np.sqrt(res.A["mid"][-n:] / np.pi) - p["rd"], "dP": res.P["in"][-n:] - res.P["out"][-n:]}
    out = {"dx": dx, "cycles": cycles, "dt": res.dt,
           "periodicity": float(np.max(np.abs(res.P["mid"][-n:] - res.P["mid"][-2 * n:-n])) / np.ptp(mine["P"])),
           "mass_error": float(res.volume_out[-1] / res.volume_in[-1] - 1)}
    d3 = np.loadtxt(base / "NumData" / "3D.txt")
    ref3 = {"P": d3[:, 1] * 1e3, "Q": d3[:, 2] * 1e-6, "dr": d3[:, 3] * 1e-3, "dP": d3[:, 4] * 1e3}
    out["vs_3D"] = all_errors(ref3, {k: on_times(d3[:, 0] + CAROTID_SHIFT_3D, t, y, period) for k, y in mine.items()})
    out["schemes_vs_3D"], out["vs_schemes"] = {}, {}
    for s in SCHEMES:
        tp, *ys = _published_cycle(base / "NumData" / f"{s}.txt", CAROTID_LAST_CYCLE, period)
        pub = dict(zip(("P", "Q", "dr", "dP"), ys))
        out["schemes_vs_3D"][s] = all_errors(ref3, {k: on_times(d3[:, 0] + CAROTID_SHIFT_3D, tp, y, period)
                                                    for k, y in pub.items()})
        out["vs_schemes"][s] = all_errors(pub, {k: on_times(tp, t, y, period) for k, y in mine.items()})
    out["mmHg"] = {"mine": [float(x / MMHG) for x in (mine["P"].min(), mine["P"].mean(), mine["P"].max())],
                   "3D": [float(x / MMHG) for x in (ref3["P"].min(), ref3["P"].mean(), ref3["P"].max())]}
    return out
