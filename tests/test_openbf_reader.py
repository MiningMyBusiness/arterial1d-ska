"""The openBF reader refuses what the solver cannot represent instead of ignoring it."""
import copy

import pytest
import yaml

import arterial1d as a1
from arterial1d.openbf import UnsupportedFeature

BASE = {
    "project_name": "t", "blood": {"rho": 1060.0, "mu": 0.004},
    "network": [
        {"label": "parent", "sn": 1, "tn": 2, "L": 0.1, "E": 4e5, "Rp": 0.01, "Rd": 0.009, "Pext": 0.0},
        {"label": "d1", "sn": 2, "tn": 3, "L": 0.1, "E": 4e5, "R0": 0.006, "R1": 1e7, "R2": 1e9, "Cc": 1e-9},
        {"label": "d2", "sn": 2, "tn": 4, "L": 0.1, "E": 4e5, "R0": 0.006, "R1": 1e7, "R2": 1e9, "Cc": 1e-9},
    ],
}


def _write(tmp_path, cfg):
    p = tmp_path / "net.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return p


def test_reads_a_small_bifurcation(tmp_path):
    net = a1.load_network(_write(tmp_path, BASE))
    assert [v.name for v in net.vessels] == ["parent", "d1", "d2"]
    assert net.n_daughters[0] == 2 and net.is_terminal.tolist() == [False, True, True]


@pytest.mark.parametrize("key,value", [("visco-elastic", True), ("h0", 0.001), ("Rt", 0.5),
                                       ("inlet_impedance_matching", True), ("Pout", 500.0), ("outlet", "wk2"), ("foo", 1)])
def test_rejects_unsupported_keys(tmp_path, key, value):
    cfg = copy.deepcopy(BASE)
    cfg["network"][1][key] = value
    with pytest.raises(UnsupportedFeature):
        a1.load_network(_write(tmp_path, cfg))


def test_rejects_anastomosis(tmp_path):
    cfg = copy.deepcopy(BASE)
    cfg["network"][2]["tn"] = 3
    with pytest.raises(UnsupportedFeature):
        a1.load_network(_write(tmp_path, cfg))


def test_rejects_outlet_without_windkessel(tmp_path):
    cfg = copy.deepcopy(BASE)
    del cfg["network"][2]["Cc"]
    with pytest.raises(UnsupportedFeature):
        a1.load_network(_write(tmp_path, cfg))


def test_rejects_two_inlets(tmp_path):
    cfg = copy.deepcopy(BASE)
    cfg["network"].append({"label": "x", "sn": 7, "tn": 8, "L": 0.1, "E": 4e5, "R0": 0.006,
                           "R1": 1e7, "R2": 1e9, "Cc": 1e-9})
    with pytest.raises(UnsupportedFeature):
        a1.load_network(_write(tmp_path, cfg))
