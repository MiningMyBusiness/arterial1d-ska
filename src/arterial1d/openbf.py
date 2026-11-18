"""Read networks in openBF's YAML format (Benemerito et al. 2024; https://github.com/INSIGNEO/openBF).

Only the subset of the format that this solver can represent is accepted. Anything else raises
UnsupportedFeature rather than being silently ignored:

  supported   label, sn, tn, L, R0 or Rp/Rd (linear taper), E, h0 (wall thickness; default the
              thickness law below), Pext (one value for the whole
              network), gamma_profile, outlets with R1, R2 and Cc (three-element Windkessel;
              `outlet: wk3` accepted);
              M, to_save, initial_pressure and initial_flow are read but not used (this solver
              sets its own mesh and initial state)
  unsupported visco-elastic walls, Rt (reflection-coefficient
              outlets), outlets without R2 and Cc, Pout != 0, inlet_impedance_matching,
              anastomoses (a node with two parents) and networks with more than one inlet

Conventions follow openBF: the wall thickness law of Boileau et al. 2015 (see
solver.wall_thickness), reference wave speed c0 = sqrt(2 E h / (3 rho r)), radii at the
reference pressure Pext, friction K_R = 2 (gamma_profile + 2) pi mu / rho.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Optional, Tuple, Union

import numpy as np
import yaml

from .solver import Network, Vessel

PathLike = Union[str, Path]

_IGNORED = {"M", "to_save", "initial_pressure", "initial_flow"}
_KNOWN = {"label", "sn", "tn", "L", "R0", "Rp", "Rd", "E", "h0", "Pext", "gamma_profile", "R1", "R2", "Cc", "Pout",
          "outlet"} | _IGNORED
_UNSUPPORTED = {"visco-elastic": "viscoelastic walls",
                "Rt": "reflection-coefficient outlets",
                "inlet_impedance_matching": "inlet impedance matching"}


class UnsupportedFeature(ValueError):
    """The network file uses part of the openBF format that this solver does not implement."""


def read_config(path: PathLike) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_network(path: PathLike, *, dx: float = 0.005, cfl: float = 0.8, junction: str = "total",
                 min_nodes: int = 3, alpha: float = 1.0) -> Network:
    """Build a Network from an openBF YAML file.

    dx: target node spacing [m]. junction: "total" (continuity of P + rho u^2/2, the default,
    which matched the published ADAN56 solution best) or "static".
    """
    if junction not in ("total", "static"):
        raise ValueError("junction must be 'total' or 'static'")
    cfg = read_config(path)
    segs = cfg["network"]
    for s in segs:
        name = s.get("label", "?")
        for k, what in _UNSUPPORTED.items():
            if k in s and s[k] not in (False, 0, 0.0):
                raise UnsupportedFeature(f"{name}: {what} ({k}) not supported")
        unknown = set(s) - _KNOWN - set(_UNSUPPORTED)
        if unknown:
            raise UnsupportedFeature(f"{name}: unknown keys {sorted(unknown)}")
        if s.get("outlet", "wk3") != "wk3":
            raise UnsupportedFeature(f"{name}: outlet type {s['outlet']!r} not supported (wk3 only)")
        if float(s.get("Pout", 0.0)) != 0.0:
            raise UnsupportedFeature(f"{name}: Pout != 0 not supported")

    by_tn: Dict[int, dict] = {}
    for s in segs:
        if s["tn"] in by_tn:
            raise UnsupportedFeature(f"node {s['tn']} has two parents (anastomosis): not supported")
        by_tn[s["tn"]] = s
    roots = [s for s in segs if s["sn"] not in by_tn]
    if len(roots) != 1:
        raise UnsupportedFeature(f"{len(roots)} inlets: exactly one is supported")
    parent_of = {s["label"]: (by_tn[s["sn"]]["label"] if s["sn"] in by_tn else None) for s in segs}

    order, placed, pending = [], set(), list(segs)          # parents before daughters, root first
    while pending:                                          # (file order otherwise)
        rest = []
        for s in pending:
            p = parent_of[s["label"]]
            if p is None or p in placed:
                order.append(s)
                placed.add(s["label"])
            else:
                rest.append(s)
        if len(rest) == len(pending):
            raise ValueError("network is not a tree")
        pending = rest

    p_ext = {float(s.get("Pext", 0.0)) for s in segs}
    if len(p_ext) != 1:
        raise UnsupportedFeature("Pext must be the same for every vessel")
    gammas = {int(s.get("gamma_profile", 2)) for s in segs}
    if len(gammas) != 1:
        raise UnsupportedFeature("gamma_profile must be the same for every vessel")

    vessels = []
    for s in order:
        if "R0" not in s and "Rp" not in s:
            raise ValueError(f"{s['label']}: missing radius")
        rp = float(s.get("Rp", s.get("R0")))
        rd = float(s.get("Rd", rp))
        v = Vessel(s["label"], float(s["L"]), rp, rd, parent_of[s["label"]], E=float(s["E"]), taper="linear",
                   h=float(s["h0"]) if "h0" in s else None)
        is_outlet = not any(parent_of[o["label"]] == s["label"] for o in segs)
        if is_outlet:
            if not all(k in s for k in ("R1", "R2", "Cc")):
                raise UnsupportedFeature(f"{s['label']}: outlets need R1, R2 and Cc (three-element Windkessel)")
            v.R1, v.R2, v.C = float(s["R1"]), float(s["R2"]), float(s["Cc"])
        vessels.append(v)
    blood = cfg["blood"]
    return Network(vessels, rho=float(blood["rho"]), mu=float(blood["mu"]), alpha=alpha,
                   friction_coeff_pi=2.0 * (gammas.pop() + 2), dx=dx, cfl=cfl, P_ven=0.0,
                   P_ref=p_ext.pop(), min_nodes=min_nodes, total_pressure=(junction == "total"))


def load_inlet(path: PathLike) -> Tuple[Callable[[np.ndarray], np.ndarray], float]:
    """An openBF inlet file (two columns: time [s], flow [m^3/s], one cycle) as a periodic
    function of time. Returns (q_in, period)."""
    d = np.loadtxt(path)
    period = float(d[-1, 0])
    return (lambda t: np.interp(np.mod(t, period), d[:, 0], d[:, 1])), period


def inlet_path(config_path: PathLike) -> Path:
    """The inlet file named in an openBF config (inlet_file, else <project_name>_inlet.dat),
    relative to the config's directory."""
    cfg = read_config(config_path)
    name = cfg.get("inlet_file", f"{cfg['project_name']}_inlet.dat")
    return Path(config_path).parent / name
