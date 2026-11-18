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


# ---- cases 2-4: single vessels and the aortic bifurcation, against 3-D (Xiao et al. 2014) -------
# Elastic tube, P = Pd + (4/3) sqrt(pi) E h / Ad (sqrt(A) - sqrt(Ad)); three-element Windkessel
# outlets; parameters from Tables III, V and VII of Boileau et al. 2015. Pd for cases 3 and 4 is
# printed as "9.4^6 kPa"; 9.466 kPa follows from the tables' Ad and A(x,0) through the tube law
# (the same calculation gives case 2's printed 10.933 kPa).
#
# phase_3d: inflow phase (s into the cycle of the inflow formula) at t = 0 of the 3-D file. Found
# by reproducing the published FVM errors (Tables IV, VI, VIII) from the files: 0.056, 0.002, 0.001.
# window: start of the last cycle in the published 1-D files (README tini). cycles: enough for a
# periodic state to 1e-7 (the Windkessel time constants of cases 3 and 4 are about 1.1 s).
# Case 3 caveat: the 3-D solution and five of the six published schemes carry a mean flow of
# 102.80-102.84 mL/s; the inflow formula in the supplement (used here, and apparently by LCG) has
# 103.085 mL/s, 0.27% more. This alone raises mean pressure by about 30 Pa (0.25% RMS).
_VESSELS = {
    "carotid": [dict(name="vessel", parent=None, L=0.126, rd=3.0e-3, h=0.3e-3, E=700.0e3,
                     R1=2.4875e8, R2=1.8697e9, C=1.7529e-10)],
    "thoracic_aorta": [dict(name="vessel", parent=None, L=0.24137, rd=1.2e-2, h=1.2e-3, E=400.0e3,
                            R1=1.1752e7, R2=1.1167e8, C=1.0163e-8)],
    "bifurcation": [dict(name="aorta", parent=None, L=0.086, rd=0.86e-2, h=1.032e-3, E=500.0e3),
                    dict(name="iliac_1", parent="aorta", L=0.085, rd=0.60e-2, h=0.72e-3, E=700.0e3,
                         R1=6.8123e7, R2=3.1013e9, C=3.6664e-10),
                    dict(name="iliac_2", parent="aorta", L=0.085, rd=0.60e-2, h=0.72e-3, E=700.0e3,
                         R1=6.8123e7, R2=3.1013e9, C=3.6664e-10)],
}
CASES = {
    "carotid": dict(folder="Benchmark2_CommonCarotidArtery", Pd=10.933e3, zeta=2, window=8.8, phase_3d=0.056, cycles=10,
                    sites={"mid": ("vessel", 0.5, "")}, dP=("vessel", 0.0, "vessel", 1.0)),
    "thoracic_aorta": dict(folder="Benchmark3_UpperThoracicAorta", Pd=9.466e3, zeta=9, window=19.045, cycles=22,
                           phase_3d=0.002, sites={"mid": ("vessel", 0.5, "")}, dP=("vessel", 0.0, "vessel", 1.0)),
    "bifurcation": dict(folder="Benchmark4_AorticBifurcation", Pd=9.466e3, zeta=9, window=28.6, phase_3d=0.001, cycles=30,
                        sites={"aorta_mid": ("aorta", 0.5, "_AortaMid"), "aorta_end": ("aorta", 1.0, "_Junction"),
                               "iliac_mid": ("iliac_1", 0.5, "_IliacMid")}, dP=None),
}
CASE_NUMBER = {"carotid": 2, "thoracic_aorta": 3, "bifurcation": 4}


def build_case(case: str, dx: float = 0.0025, junction: str = "total") -> Network:
    c = CASES[case]
    vs = [Vessel(v["name"], v["L"], v["rd"], v["rd"], v["parent"], E=v["E"], h=v["h"],
                 R1=v.get("R1"), R2=v.get("R2"), C=v.get("C")) for v in _VESSELS[case]]
    return Network(vs, rho=1060.0, mu=4.0e-3, alpha=1.0, friction_coeff_pi=2.0 * (c["zeta"] + 2), dx=dx,
                   P_ven=0.0, P_ref=c["Pd"], total_pressure=(junction == "total"))


def _published_cycle(path, t0, period):
    """One cycle of a published solution: (t in cycle, {P [Pa], Q [m^3/s], dr [m], dP [Pa]})."""
    d = np.loadtxt(path)
    d = d[(d[:, 0] >= t0 - 1e-9) & (d[:, 0] < t0 + period - 1e-9)]
    scale = (1e3, 1e-6, 1e-3, 1e3)
    return d[:, 0] - t0, {k: d[:, j + 1] * scale[j] for j, k in enumerate(("P", "Q", "dr", "dP")[:d.shape[1] - 1])}


def run_case(case: str, root: Optional[os.PathLike] = None, *, dx: float = 0.0025, cycles: Optional[int] = None,
             junction: str = "total") -> dict:
    """Run case 2, 3 or 4 and compare, at each published site, P, Q, dr (and the inlet-outlet
    pressure drop for single vessels) with the 3-D solution and with the six published 1-D
    schemes, using the paper's error measures (%)."""
    r = reference_root(root)
    if r is None:
        raise FileNotFoundError(f"Boileau et al. 2015 supplementary data not found; set {ENV}")
    c = CASES[case]
    cycles = cycles or c["cycles"]
    base = r / c["folder"]
    f, period, _ = read_inflow_formula(base / "inflow.txt")
    net = build_case(case, dx, junction)
    probes = {k: net.node(v, x) for k, (v, x, _) in c["sites"].items()}
    if c["dP"]:
        probes["dP_in"], probes["dP_out"] = net.node(c["dP"][0], c["dP"][1]), net.node(c["dP"][2], c["dP"][3])
    res = net.simulate(lambda t: f(t) * 1e-6, period, cycles, probes, P_init=c["Pd"])
    if not res.finite:
        raise RuntimeError("simulation did not stay finite")
    n = int(round(period * 1000))
    t = np.mod(res.t[-n:] + 1e-9, period) - 1e-9          # inflow phase (true time mod period)
    vessel_of = {k: v for k, (v, _, _) in c["sites"].items()}
    rd = {v["name"]: v["rd"] for v in _VESSELS[case]}
    out = {"case": CASE_NUMBER[case], "dx": dx, "cycles": cycles, "dt": res.dt, "junction": junction,
           "periodicity": float(max(np.max(np.abs(res.P[k][-n:] - res.P[k][-2 * n:-n])) / np.ptp(res.P[k][-n:])
                                    for k in c["sites"])),
           "mass_error": float(res.volume_out[-1] / res.volume_in[-1] - 1), "sites": {}}
    for site, (vname, _, suffix) in c["sites"].items():
        mine = {"P": res.P[site][-n:], "Q": res.Q[site][-n:], "dr": np.sqrt(res.A[site][-n:] / np.pi) - rd[vname]}
        if c["dP"]:
            mine["dP"] = res.P["dP_in"][-n:] - res.P["dP_out"][-n:]
        d3 = np.loadtxt(base / "NumData" / f"3D{suffix}.txt")
        t3 = d3[:, 0] + c["phase_3d"]
        _, ref3 = _published_cycle(base / "NumData" / f"3D{suffix}.txt", 0.0, period)
        so = {"vs_3D": all_errors(ref3, {k: on_times(t3, t, y, period) for k, y in mine.items()}),
              "schemes_vs_3D": {}, "vs_schemes": {},
              "mmHg": {"mine": [float(x / MMHG) for x in (mine["P"].min(), mine["P"].mean(), mine["P"].max())],
                       "3D": [float(x / MMHG) for x in (ref3["P"].min(), ref3["P"].mean(), ref3["P"].max())]}}
        win_phase = np.mod(c["window"], period)
        for s in SCHEMES:
            tp, pub = _published_cycle(base / "NumData" / f"{s}{suffix}.txt", c["window"], period)
            if not all(np.isfinite(y).all() for y in pub.values()):
                continue                                   # e.g. case 4 LCG files contain NaN
            tp = tp + win_phase                            # published window start -> inflow phase
            so["schemes_vs_3D"][s] = all_errors(ref3, {k: on_times(t3, tp, y, period) for k, y in pub.items()})
            so["vs_schemes"][s] = all_errors(pub, {k: on_times(tp, t, y, period) for k, y in mine.items() if k in pub})
        out["sites"][site] = so
    return out


def carotid(root: Optional[os.PathLike] = None, *, dx: float = 0.0025, cycles: int = 10) -> dict:
    """Case 2 (kept for the original interface): run_case('carotid') with its single site flattened."""
    o = run_case("carotid", root, dx=dx, cycles=cycles)
    site = o.pop("sites")["mid"]
    o.update(site)
    return o


# ---- case 5: 37-artery in-vitro network (Matthys et al. 2007; Boileau et al. 2015, section 4.5) ----
# Geometry, wall thickness and terminal resistances from the supplement's mesh.txt (resistances
# there are in units of 1e9 Pa s/m^3); rho 1050, mu 2.5 mPa s, zeta 9, E 1.2 MPa, Pd = Pext = 0
# (A0 at zero pressure), outflow pressure 432.6 Pa, linear radius taper (Table IX).
# Measurement sites (not listed in the paper) identified by matching every segment and position
# against the published FVM output: all are segment midpoints, each matching to 0.2-0.5% with the
# next candidate at least 4x worse.
NETWORK37_SITES = {"AorticArchII": 10, "ThoracicAortaII": 17, "LSubclavianI": 11, "RIliacFemoralII": 29,
                   "LUlnar": 14, "RAnteriorTibial": 34, "RUlnar": 7, "Splenic": 20}
NETWORK37_PUBLISHED_END = 15.0     # the published 1-D runs end at 14.999 s
NETWORK37_OFFSET = -1e-3           # 1 ms time-label offset; with it the code reproduces Table X (FVM)


def build_network37(root: Optional[os.PathLike] = None, dx: float = 0.0025, junction: str = "total") -> Network:
    r = reference_root(root)
    if r is None:
        raise FileNotFoundError(f"Boileau et al. 2015 supplementary data not found; set {ENV}")
    m = np.loadtxt(r / "Benchmark5_37ArteryNetwork" / "mesh.txt")
    seg_ending_at = {int(row[2]): int(row[0]) for row in m}
    starts = set(m[:, 1].astype(int))
    vs = []
    for row in m:
        k, n0, n1 = int(row[0]), int(row[1]), int(row[2])
        parent = seg_ending_at.get(n0)
        terminal = n1 not in starts
        vs.append(Vessel(f"s{k}", row[3], row[4], row[5], None if parent is None else f"s{parent}", E=1.2e6,
                         h=row[6], taper="linear", R1=row[9] * 1e9 if terminal else None))
    return Network(vs, rho=1050.0, mu=2.5e-3, alpha=1.0, friction_coeff_pi=22.0, dx=dx, P_ven=432.6, P_ref=0.0,
                   total_pressure=(junction == "total"))


def _time_windows(path):
    import re
    raw = Path(path).read_bytes().decode("utf-16")
    return {a: float(b) for a, b in re.findall(r"(\w+)\s*tini\s*=\s*([0-9.]+)", raw)}


def network37(root: Optional[os.PathLike] = None, *, dx: float = 0.0025, cycles: int = 30,
              junction: str = "total") -> dict:
    """Run case 5 and compare P and Q at the eight measurement sites with the in-vitro data and
    with the six published schemes (the paper's error measures, %)."""
    r = reference_root(root)
    if r is None:
        raise FileNotFoundError(f"Boileau et al. 2015 supplementary data not found; set {ENV}")
    base = r / "Benchmark5_37ArteryNetwork"
    f, period, _ = read_inflow_formula(base / "inflow.txt")
    net = build_network37(r, dx, junction)
    res = net.simulate(lambda t: f(t) * 1e-6, period, cycles,
                       {nm: net.node(f"s{k}", 0.5) for nm, k in NETWORK37_SITES.items()}, P_init=0.0)
    if not res.finite:
        raise RuntimeError("simulation did not stay finite")
    n = int(round(period * 1000))
    t = np.mod(res.t[-n:] + 1e-9, period) - 1e-9
    tini = _time_windows(base / "time_shifts.txt")
    pub = {s: (np.loadtxt(base / "NumData" / f"{s}_P.txt"), np.loadtxt(base / "NumData" / f"{s}_Q.txt"))
           for s in SCHEMES}
    out = {"case": 5, "dx": dx, "cycles": cycles, "dt": res.dt, "junction": junction,
           "periodicity": float(max(np.max(np.abs(res.P[k][-n:] - res.P[k][-2 * n:-n])) / np.ptp(res.P[k][-n:])
                                    for k in NETWORK37_SITES)),
           "mass_error": float(res.volume_out[-1] / res.volume_in[-1] - 1), "sites": {}}
    for j, nm in enumerate(NETWORK37_SITES):
        ep, eq = np.loadtxt(base / "ExpData" / f"{nm}_P.txt"), np.loadtxt(base / "ExpData" / f"{nm}_Q.txt")
        te = ep[:, 0] + tini[nm] + NETWORK37_OFFSET
        keep = te < NETWORK37_PUBLISHED_END - 1e-9        # as the paper: samples inside the published run
        te = te[keep]
        ref = {"P": ep[keep, 1] * 1e3, "Q": eq[keep, 1] * 1e-6}
        mine = {"P": on_times(te, t, res.P[nm][-n:], period), "Q": on_times(te, t, res.Q[nm][-n:], period)}
        so = {"n_samples": int(keep.sum()), "vs_exp": all_errors(ref, mine), "schemes_vs_exp": {}, "vs_schemes": {},
              "mmHg": {"mine": [float(x / MMHG) for x in (mine["P"].min(), mine["P"].mean(), mine["P"].max())],
                       "exp": [float(x / MMHG) for x in (ref["P"].min(), ref["P"].mean(), ref["P"].max())]}}
        for s, (P, Q) in pub.items():
            col = j + 1
            if not (np.isfinite(P[:, col]).all() and np.isfinite(Q[:, col]).all()):
                continue
            so["schemes_vs_exp"][s] = all_errors(ref, {"P": np.interp(te, P[:, 0], P[:, col] * 1e3),
                                                       "Q": np.interp(te, Q[:, 0], Q[:, col] * 1e-6)})
            tp = np.mod(P[-n:, 0] + 1e-9, period) - 1e-9   # last cycle of the published run
            so["vs_schemes"][s] = all_errors({"P": P[-n:, col] * 1e3, "Q": Q[-n:, col] * 1e-6},
                                             {"P": on_times(tp, t, res.P[nm][-n:], period),
                                              "Q": on_times(tp, t, res.Q[nm][-n:], period)})
        out["sites"][nm] = so
    return out
