#!/usr/bin/env python
"""Compare arterial1d with the published ADAN56 solutions of Boileau et al. 2015 (benchmark 6).

    python scripts/benchmark_adan56.py --data /path/to/SuppMaterial_1DBenchmark_Boileau_etal
    python scripts/benchmark_adan56.py --reference FEM --junction static --json out.json

Prints, per site, the pressure RMS difference as % of pulse pressure and the flow RMS difference
as % of peak flow, with (min/mean/max) pressures of both solutions.
"""
import argparse
import json

from arterial1d import benchmarks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", help=f"supplementary data folder (default: ${benchmarks.ENV})")
    ap.add_argument("--dx", type=float, default=0.0025)
    ap.add_argument("--cycles", type=int, default=15)
    ap.add_argument("--junction", choices=["total", "static"], default="total")
    ap.add_argument("--reference", choices=benchmarks.SCHEMES, default="FVM")
    ap.add_argument("--json", help="also write the results to this file")
    a = ap.parse_args()
    res = benchmarks.adan56(a.data, dx=a.dx, cycles=a.cycles, junction=a.junction, reference=a.reference)
    print(f"ADAN56, dx = {a.dx} m, {a.cycles} cycles, {a.junction}-pressure junctions, against {a.reference}\n")
    print(f"{'site':26s} {'P, mmHg (min/mean/max)':>24s} {'reference':>20s} {'P RMS %PP':>10s} {'Q RMS %peak':>12s}")
    for s, v in res.items():
        p, r = v["P_mmHg"], v["P_ref_mmHg"]
        print(f"{s:26s} {p[0]:7.1f}/{p[1]:6.1f}/{p[2]:6.1f}   {r[0]:6.1f}/{r[1]:6.1f}/{r[2]:6.1f} "
              f"{100 * v['P_rms_over_PP']:10.2f} {100 * v['Q_rms_over_peak']:12.2f}")
    ps = sorted(v["P_rms_over_PP"] for v in res.values())
    print(f"\npressure RMS, % of pulse pressure: median {100 * ps[len(ps) // 2]:.2f}, worst {100 * ps[-1]:.2f}")
    if a.json:
        with open(a.json, "w") as f:
            json.dump({"args": vars(a), "sites": res}, f, indent=2)


if __name__ == "__main__":
    main()
