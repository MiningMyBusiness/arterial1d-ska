"""Regression against the published ADAN56 solution (Boileau et al. 2015, benchmark 6).

Skipped unless the supplementary data are available (set ARTERIAL1D_BOILEAU_DIR; see README).
"""
import numpy as np
import pytest

from arterial1d import benchmarks

ROOT = benchmarks.reference_root()
needs_data = pytest.mark.skipif(ROOT is None, reason=f"set {benchmarks.ENV} to the Boileau et al. 2015 data")


@needs_data
def test_adan56_matches_published_finite_volume_solution():
    """Within 2% of pulse pressure at all 12 sites (the published schemes differ from FVM by
    0.2-4.5% at their worst site). At dx = 2.5 mm this solver gives 0.5% median, 1.1% worst."""
    res = benchmarks.adan56(ROOT)
    worst = max(v["P_rms_over_PP"] for v in res.values())
    assert worst < 0.02, {a: round(100 * v["P_rms_over_PP"], 2) for a, v in res.items()}
    for v in res.values():
        assert np.allclose(v["P_mmHg"], v["P_ref_mmHg"], atol=1.0)


TABLE_IV_FVM = {"P": (0.28, 0.66, -0.29, 0.29), "Q": (0.29, 1.07, -0.58, 0.24),
                "dr": (1.07, 2.28, -1.68, 0.12), "dP": (4.21, 15.60, -14.83, 4.91)}


@needs_data
def test_metric_code_reproduces_table_iv_for_published_fvm():
    """The error measures and the 3-D time alignment reproduce Boileau et al. 2015 Table IV
    (FVM column) from the published files alone: all 16 entries within 0.015 percentage points
    (13 round to the printed value; DIAS Q 0.249 vs 0.24, DIAS dr 0.111 vs 0.12, MAX dP 15.611 vs 15.60)."""
    o = benchmarks.carotid(ROOT, dx=0.005, cycles=2)
    for q, vals in TABLE_IV_FVM.items():
        for m, v in zip(("RMS", "MAX", "SYS", "DIAS"), vals):
            assert o["schemes_vs_3D"]["FVM"][q][m] == pytest.approx(v, abs=0.015), (q, m)


@needs_data
def test_carotid_agrees_with_3d_as_well_as_the_published_schemes():
    """Case 2. Published schemes vs 3-D (Table IV), RMS %: P 0.20-0.28, Q 0.26-0.39, dr 0.96-1.07,
    dP 4.21-4.47. arterial1d at dx = 2.5 mm: 0.24, 0.28, 1.00, 4.27."""
    o = benchmarks.carotid(ROOT)
    v = o["vs_3D"]
    assert v["P"]["RMS"] < 0.30 and v["Q"]["RMS"] < 0.40 and v["dr"]["RMS"] < 1.10 and v["dP"]["RMS"] < 4.5
    assert o["vs_schemes"]["FVM"]["P"]["RMS"] < 0.2 and o["vs_schemes"]["FVM"]["Q"]["RMS"] < 0.2
    assert o["periodicity"] < 1e-6 and abs(o["mass_error"]) < 1e-5


# Published FVM errors vs 3-D (RMS P, MAX P, RMS Q, MAX Q, RMS dr, MAX dr): Tables VI and VIII
TABLE_VI_VIII_FVM = {
    ("thoracic_aorta", "mid"): (1.14, 3.18, 2.17, 7.07, 2.44, 7.20),
    ("bifurcation", "aorta_mid"): (0.40, 0.68, 0.81, 2.44, 2.49, 4.01),
    ("bifurcation", "aorta_end"): (0.44, 0.80, 1.08, 3.28, 4.08, 6.95),
    ("bifurcation", "iliac_mid"): (0.47, 0.92, 0.68, 1.80, 4.31, 7.48),
}


@pytest.fixture(scope="module")
def cases_3_4():
    return {c: benchmarks.run_case(c, ROOT) for c in ("thoracic_aorta", "bifurcation")} if ROOT else None


@needs_data
@pytest.mark.parametrize("case,site", list(TABLE_VI_VIII_FVM))
def test_metric_code_reproduces_tables_vi_viii_for_published_fvm(cases_3_4, case, site):
    e = cases_3_4[case]["sites"][site]["schemes_vs_3D"]["FVM"]
    got = (e["P"]["RMS"], e["P"]["MAX"], e["Q"]["RMS"], e["Q"]["MAX"], e["dr"]["RMS"], e["dr"]["MAX"])
    assert got == pytest.approx(TABLE_VI_VIII_FVM[(case, site)], abs=0.035)


@needs_data
def test_thoracic_aorta_against_3d(cases_3_4):
    """Case 3. Published schemes vs 3-D, RMS %: P 1.08-1.21, Q 2.17-2.55, dr 2.32-2.90, dP 7.14-7.70.
    arterial1d: 1.22, 2.22, 2.43, 7.23 (pressure slightly above the range because the supplement's
    inflow formula carries 0.27% more mean flow than the 3-D solution; see benchmarks.py)."""
    o = cases_3_4["thoracic_aorta"]
    v = o["sites"]["mid"]["vs_3D"]
    assert v["P"]["RMS"] < 1.30 and v["Q"]["RMS"] < 2.6 and v["dr"]["RMS"] < 2.95 and v["dP"]["RMS"] < 7.8
    assert o["periodicity"] < 1e-6


@needs_data
def test_bifurcation_against_3d_and_published_schemes(cases_3_4):
    """Case 4. arterial1d matches the published McC and FEM solutions to about 0.01% RMS pressure
    at all three sites, and is within the published range against 3-D."""
    o = cases_3_4["bifurcation"]
    upper = {"aorta_mid": (0.45, 0.95, 2.6), "aorta_end": (0.50, 1.25, 4.2), "iliac_mid": (0.52, 0.95, 4.4)}
    for site, (p, q, dr) in upper.items():
        so = o["sites"][site]
        v = so["vs_3D"]
        assert v["P"]["RMS"] < p and v["Q"]["RMS"] < q and v["dr"]["RMS"] < dr, site
        assert so["vs_schemes"]["McC"]["P"]["RMS"] < 0.05 and so["vs_schemes"]["FEM"]["P"]["RMS"] < 0.05, site
    assert o["periodicity"] < 1e-6
