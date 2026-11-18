"""Error measures of Boileau et al. 2015 (Eqs. 10-13), in percent, over one cardiac cycle.

Pressure errors are relative to the reference value at each instant; errors in flow Q, change in
radius dr and pressure drop dP are relative to the maximum of the reference over the cycle (the
paper: flow is normalised by its maximum, and dr and dP use "the same metrics as for the flow
rate"). The 1-D solution is interpolated, periodically, onto the reference's sample times.
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def on_times(t_ref, t, y, period):
    o = np.argsort(np.mod(t, period))
    return np.interp(np.mod(t_ref, period), np.mod(t, period)[o], np.asarray(y)[o], period=period)


def errors(ref: np.ndarray, sol: np.ndarray, relative_pointwise: bool) -> Dict[str, float]:
    """ref, sol: the same quantity at the same instants (one cycle)."""
    if relative_pointwise:
        e = (sol - ref) / ref
        dias = (sol.min() - ref.min()) / ref.min()
    else:
        scale = np.max(ref)
        e = (sol - ref) / scale
        dias = (sol.min() - ref.min()) / scale
    sys_ = (sol.max() - ref.max()) / ref.max()
    return {"RMS": 100 * float(np.sqrt(np.mean(e ** 2))), "MAX": 100 * float(np.max(np.abs(e))),
            "SYS": 100 * float(sys_), "DIAS": 100 * float(dias)}


def all_errors(ref: Dict[str, np.ndarray], sol: Dict[str, np.ndarray]) -> Dict[str, Dict[str, float]]:
    """ref/sol keyed by 'P', 'Q', 'dr', 'dP' (any subset present in both)."""
    return {k: errors(ref[k], sol[k], relative_pointwise=(k == "P")) for k in ("P", "Q", "dr", "dP")
            if k in ref and k in sol}
