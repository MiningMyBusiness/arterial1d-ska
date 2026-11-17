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
from .solver import MMHG

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
