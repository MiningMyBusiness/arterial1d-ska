#!/usr/bin/env python
"""Boileau et al. 2015 cases 2-4 against the 3-D solutions (ADAN56: benchmark_adan56.py).

    python scripts/benchmark_boileau.py carotid        --data /path/to/SuppMaterial_1DBenchmark_Boileau_etal
    python scripts/benchmark_boileau.py thoracic_aorta
    python scripts/benchmark_boileau.py bifurcation
    python scripts/benchmark_boileau.py network37      # case 5, against the in-vitro measurements

For each published site: the paper's error measures (Eqs. 10-13, %) for arterial1d against 3-D,
next to the published 1-D schemes recomputed the same way, then arterial1d against each scheme.
"""
import argparse
import json

from arterial1d import benchmarks

ROWS = [(q, m) for m in ("RMS", "MAX", "SYS", "DIAS") for q in ("P", "Q", "dr", "dP")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=list(benchmarks.CASES) + ["network37"])
    ap.add_argument("--data", help=f"supplementary data folder (default: ${benchmarks.ENV})")
    ap.add_argument("--dx", type=float, default=0.0025)
    ap.add_argument("--cycles", type=int, help="default: per case")
    ap.add_argument("--junction", choices=["total", "static"], default="total")
    ap.add_argument("--json")
    a = ap.parse_args()
    if a.case == "network37":
        o = benchmarks.network37(a.data, dx=a.dx, cycles=a.cycles or 30, junction=a.junction)
        ref, key = "in vitro", "vs_exp"
    else:
        o = benchmarks.run_case(a.case, a.data, dx=a.dx, cycles=a.cycles, junction=a.junction)
        ref, key = "3-D", "vs_3D"
    print(f"case {o['case']} ({a.case}): dx = {a.dx} m, dt = {o['dt'] * 1e6:.1f} us, {o['cycles']} cycles, "
          f"periodicity {o['periodicity']:.1e}, mass error {o['mass_error']:+.1e}")
    for site, so in o["sites"].items():
        p = so["mmHg"]
        r = p.get("3D", p.get("exp"))
        print(f"\n== {site}: pressure min/mean/max arterial1d {p['mine'][0]:.1f}/{p['mine'][1]:.1f}/{p['mine'][2]:.1f}, "
              f"{ref} {r[0]:.1f}/{r[1]:.1f}/{r[2]:.1f} mmHg")
        names = list(so["schemes_" + key])
        print(f"errors against the {ref} data, %   (published schemes recomputed with the same code)")
        print(f"{'':10s} {'arterial1d':>10s}" + "".join(f"{s:>8s}" for s in names))
        for q, m in ROWS:
            if q in so[key]:
                print(f"{m + ' ' + q:10s} {so[key][q][m]:10.2f}" +
                      "".join(f"{so['schemes_' + key][s][q][m]:8.2f}" for s in names))
        print("arterial1d against each published scheme, RMS %")
        for q in ("P", "Q", "dr", "dP"):
            if q in so[key]:
                print(f"{'RMS ' + q:10s} {'':10s}" + "".join(f"{so['vs_schemes'][s][q]['RMS']:8.2f}" for s in names))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(o, f, indent=2)


if __name__ == "__main__":
    main()
