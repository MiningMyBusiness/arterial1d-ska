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
