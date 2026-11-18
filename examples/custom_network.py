#!/usr/bin/env python
"""Build a network in code: an aorta that bifurcates into two iliac arteries, Windkessel outlets.

    python examples/custom_network.py

Vessels are listed root first, each after its parent. Stiffness from a Young's modulus E and a wall
thickness h (or the empirical thickness law if h is omitted, or wave speeds c_in/c_out). Radii are
at the reference pressure P_ref. Outlets: three-element Windkessel (R1, R2, C), or a pure
resistance R1 (R2 = C = None); R1 = None matches the characteristic impedance.
"""
import numpy as np

import arterial1d as a1
from arterial1d import MMHG, Network, Vessel

vessels = [
    Vessel("aorta", 0.086, 0.0086, 0.0086, None, E=500e3, h=1.032e-3),
    Vessel("iliac_L", 0.085, 0.0060, 0.0060, "aorta", E=700e3, h=0.72e-3, R1=6.8123e7, R2=3.1013e9, C=3.6664e-10),
    Vessel("iliac_R", 0.085, 0.0060, 0.0060, "aorta", E=700e3, h=0.72e-3, R1=6.8123e7, R2=3.1013e9, C=3.6664e-10),
]
net = Network(vessels, rho=1060.0, mu=4e-3, friction_coeff_pi=22.0,   # 22 pi: flat profile, zeta = 9
              P_ref=9.466e3)                                          # radii above are at 9.466 kPa

period = 1.1
def q_in(t):                                                          # a smooth systolic pulse, 8 mL/s mean
    tc = np.mod(t, period)
    return np.where(tc < 0.35, 8e-6 * period / 0.35 * np.pi / 2 * np.sin(np.pi * tc / 0.35), 0.0)

res = net.simulate(q_in, period, n_cycles=30, P_init=net.P_ref,
                   probes={"aorta mid": net.node("aorta", 0.5), "iliac mid": net.node("iliac_L", 0.5)})
n = int(round(period * 1000))
print(f"dt {res.dt * 1e6:.0f} us, max Courant {res.max_courant:.2f}, periodic to "
      f"{np.max(np.abs(res.P['aorta mid'][-n:] - res.P['aorta mid'][-2 * n:-n])) / np.ptp(res.P['aorta mid'][-n:]):.1e}")
for k in ("aorta mid", "iliac mid"):
    p = res.P[k][-n:] / MMHG
    print(f"{k:10s} pressure {p.min():.1f} / {p.mean():.1f} / {p.max():.1f} mmHg, "
          f"mean flow {res.Q[k][-n:].mean() * 1e6:.2f} mL/s")
