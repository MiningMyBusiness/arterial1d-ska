#!/usr/bin/env python
"""Quickstart: the bundled ADAN56 network (56 arteries, 77 segments) with its benchmark inflow.

    python examples/quickstart_adan56.py            # writes quickstart_adan56.png if matplotlib is installed

Any openBF network file within the supported subset loads the same way:
    net = a1.load_network("my_network.yaml")
    q_in, period = a1.load_inlet(a1.openbf.inlet_path("my_network.yaml"))
"""
import numpy as np

import arterial1d as a1

net = a1.load_network(a1.adan56_path())                     # dx = 2.5 mm, total-pressure junctions
q_in, period = a1.load_inlet(a1.openbf.inlet_path(a1.adan56_path()))
sites = {"aortic arch": ("aortic_arch_I", 0.5), "radial (wrist)": ("radial_L", 0.95),
         "femoral": ("femoral_R_II", 0.5), "anterior tibial": ("anterior_tibial_R", 0.5)}
probes = {k: net.node(v, x) for k, (v, x) in sites.items()}

res = net.simulate(q_in, period, n_cycles=10, probes=probes, P_init=net.P_ref)
n = int(round(period * 1000))                               # samples per cycle (recorded at 1 kHz)
print(f"time step {res.dt * 1e6:.0f} us, max Courant number {res.max_courant:.2f}, "
      f"mass balance {100 * (res.volume_out[-1] / res.volume_in[-1] - 1):+.3f}%")
for k in sites:
    p = res.P[k][-n:] / a1.MMHG
    print(f"{k:16s} {p.min():6.1f} / {p.mean():6.1f} / {p.max():6.1f} mmHg (min / mean / max)")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("(install matplotlib for the figure)")
fig, ax = plt.subplots(figsize=(6.4, 3.6))
t = np.arange(n) / 1000.0
for k, color in zip(sites, ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")):
    ax.plot(t, res.P[k][-n:] / a1.MMHG, color=color, lw=2, label=k)
ax.set_xlabel("time in cycle (s)")
ax.set_ylabel("pressure (mmHg)")
ax.grid(color="#e6e5e0", lw=0.8)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("quickstart_adan56.png", dpi=150)
print("wrote quickstart_adan56.png")
