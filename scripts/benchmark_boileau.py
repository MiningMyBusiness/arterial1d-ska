#!/usr/bin/env python
"""Boileau et al. 2015 benchmark cases other than ADAN56 (see benchmark_adan56.py).

    python scripts/benchmark_boileau.py carotid --data /path/to/SuppMaterial_1DBenchmark_Boileau_etal

Prints the paper's error measures (Eqs. 10-13, %) for arterial1d against the 3-D solution, next to
the published 1-D schemes computed the same way, and arterial1d against each published scheme.
"""
import argparse
import json

from arterial1d import benchmarks

ROWS = [(q, m) for m in ("RMS", "MAX", "SYS", "DIAS") for q in ("P", "Q", "dr", "dP")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=["carotid"])
    ap.add_argument("--data", help=f"supplementary data folder (default: ${benchmarks.ENV})")
    ap.add_argument("--dx", type=float, default=0.0025)
    ap.add_argument("--cycles", type=int, default=10)
    ap.add_argument("--json")
    a = ap.parse_args()
    o = benchmarks.carotid(a.data, dx=a.dx, cycles=a.cycles)
    print(f"case 2, common carotid: dx = {a.dx} m, dt = {o['dt'] * 1e6:.1f} us, {a.cycles} cycles, "
          f"periodicity {o['periodicity']:.1e}, mass error {o['mass_error']:+.1e}")
    p = o["mmHg"]
    print(f"mid-vessel pressure min/mean/max: arterial1d {p['mine'][0]:.1f}/{p['mine'][1]:.1f}/{p['mine'][2]:.1f}, "
          f"3-D {p['3D'][0]:.1f}/{p['3D'][1]:.1f}/{p['3D'][2]:.1f} mmHg\n")
    names = list(o["schemes_vs_3D"])
    print("errors against the 3-D solution, %   (published schemes recomputed with the same code)")
    print(f"{'':10s} {'arterial1d':>10s}" + "".join(f"{s:>8s}" for s in names))
    for q, m in ROWS:
        print(f"{m + ' ' + q:10s} {o['vs_3D'][q][m]:10.2f}" + "".join(f"{o['schemes_vs_3D'][s][q][m]:8.2f}" for s in names))
    print("\narterial1d against each published scheme, RMS %")
    print(f"{'':10s}" + "".join(f"{s:>8s}" for s in names))
    for q in ("P", "Q", "dr", "dP"):
        print(f"{'RMS ' + q:10s}" + "".join(f"{o['vs_schemes'][s][q]['RMS']:8.2f}" for s in names))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(o, f, indent=2)


if __name__ == "__main__":
    main()
