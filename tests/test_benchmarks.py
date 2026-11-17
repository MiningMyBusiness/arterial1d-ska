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
