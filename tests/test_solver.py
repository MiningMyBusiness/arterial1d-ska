"""Physics checks of the solver on a uniform tube and on ADAN56."""
import numpy as np
import pytest

import arterial1d as a1
from arterial1d import MMHG, Network, Vessel

ADAN56 = a1.adan56_path()


def _uniform_tube(c=6.0, r=0.01, length=1.0, dx=0.0025):
    # huge C: the Windkessel pressure stays put, so the outlet is a pure matched resistance
    return Network([Vessel("tube", length, r, r, None, c_in=c, c_out=c, R1=None, R2=1e12, C=1.0)], mu=0.0, dx=dx)


def _pulse(t):
    return 1e-6 * np.exp(-0.5 * ((t - 0.03) / 0.006) ** 2)        # small pulse: linear regime


def test_pulse_travels_at_the_tube_law_wave_speed():
    c = 6.0
    net = _uniform_tube(c=c)
    r = net.simulate(_pulse, 0.4, 1, {"a": net.node("tube", 0.2), "b": net.node("tube", 0.8)}, P_init=0.0)
    t_a, t_b = r.t[np.argmax(r.Q["a"])], r.t[np.argmax(r.Q["b"])]
    assert 0.6 / (t_b - t_a) == pytest.approx(c, rel=0.02)


def test_matched_outlet_reflects_little():
    net = _uniform_tube()
    r = net.simulate(_pulse, 0.4, 1, {"mid": net.node("tube", 0.5)}, P_init=0.0)
    p = r.P["mid"]
    incident = p[r.t < 0.17].max()
    reflected = np.abs(p[r.t > 0.2]).max()          # back at mid-tube after 0.5 m out and 0.5 m back
    assert reflected < 0.03 * incident


@pytest.fixture(scope="module")
def adan56():
    return a1.load_network(ADAN56)


def test_rest_state_is_preserved_on_the_tapered_network(adan56):
    net = adan56
    pr = {"root": 0, "radial": net.node("radial_L", 0.95), "tibial": net.node("anterior_tibial_R", 0.5)}
    r = net.simulate(lambda t: np.zeros_like(t), 0.5, 2, pr, P_init=net.P_ven)
    assert r.finite
    assert max(np.abs(q).max() for q in r.Q.values()) < 1e-12
    assert all(np.allclose(p, net.P_ven, atol=1e-6) for p in r.P.values())


def test_grid_spans_vessel_length(adan56):
    net = adan56
    for v, x in zip(net.vessels, net.x):
        assert x[0] == 0.0 and x[-1] == pytest.approx(v.length)
    i = net.index["radial_L"]
    assert np.sqrt(net.A0[net.end[i] - 1] / np.pi) == pytest.approx(net.vessels[i].r_out)


def test_network_is_adan56_as_published(adan56):
    src = a1.openbf.read_config(ADAN56)["network"]
    net = adan56
    assert len(net.vessels) == 77 and int(net.is_terminal.sum()) == 31
    r_tot = 1.0 / sum(1.0 / (float(s["R1"]) + float(s["R2"])) for s in src if "R1" in s)
    mine = 1.0 / np.sum(1.0 / (net.R1 + net.R2)[net.is_terminal])
    assert mine == pytest.approx(r_tot, rel=1e-12)
    assert net.K_R == pytest.approx(8 * np.pi * 0.004 / 1060.0)    # gamma_profile 2
    c0 = net.wave_speed0()
    assert 3.5 < c0[0] < 6.0                      # aortic root wave speed, m/s
    assert c0[net.node("radial_L", 1.0)] > 1.5 * c0[0]   # stiffer towards the periphery


@pytest.fixture(scope="module")
def benchmark_run(adan56):
    """ADAN56 with its benchmark inflow (Boileau et al. 2015, case 6), 10 cycles."""
    net = adan56
    q, period = a1.load_inlet(a1.openbf.inlet_path(ADAN56))
    pr = {"central": net.node("aortic_arch_I", 0.5), "wrist": net.node("radial_L", 0.95),
          "femoral": net.node("femoral_L_II", 0.5)}
    return net.simulate(q, period, 10, pr, P_init=10000.0), period


def test_mass_is_conserved(benchmark_run):
    """Volume in equals volume out per cycle within 0.5% at the periodic state."""
    r, _ = benchmark_run
    assert r.finite
    assert r.volume_out[-1] == pytest.approx(r.volume_in[-1], rel=0.005)


def test_periodic_steady_state(benchmark_run):
    r, period = benchmark_run
    n = int(round(period * 1000))
    for key in ("central", "wrist"):
        last, prev = r.P[key][-n:], r.P[key][-2 * n:-n]
        assert np.max(np.abs(last - prev)) < 0.005 * (last.max() - last.min())


def test_pressures_are_physiological(benchmark_run):
    r, period = benchmark_run
    n = int(round(period * 1000))
    p = r.P["central"][-n:] / MMHG
    assert 70 < p.min() < 95 and 110 < p.max() < 135 and 95 < p.mean() < 115
    wrist = r.P["wrist"][-n:] / MMHG
    assert 0 < p.mean() - wrist.mean() < 25                 # mean pressure falls along the arm
    fem = r.P["femoral"][-n:] / MMHG
    assert fem.max() - fem.min() > 1.3 * (p.max() - p.min())   # pulse pressure amplifies towards the leg


def test_prescribed_wall_thickness_sets_wave_speed():
    E, h, r, rho = 7e5, 0.3e-3, 3e-3, 1060.0
    net = Network([Vessel("v", 0.1, r, r, None, E=E, h=h, R1=1e8, R2=1e9, C=1e-10)], rho=rho)
    assert net.wave_speed0()[0] == pytest.approx(np.sqrt(2 * E * h / (3 * rho * r)), rel=1e-12)


def test_resistance_outlet_gives_ohmic_pressure_drop_in_steady_flow():
    """Steady 5 mL/s through a short viscous tube into R = 1e8 Pa s/m^3 with P_ven = 432.6 Pa:
    outlet pressure P_ven + R Q, and the Windkessel state does not move."""
    R, Pv, q = 1e8, 432.6, 5e-6
    net = Network([Vessel("v", 0.2, 3e-3, 3e-3, None, c_in=6.0, c_out=6.0, R1=R)], mu=0.0035, P_ven=Pv, dx=0.005)
    i = net.node("v", 1.0)
    r = net.simulate(lambda t: np.full_like(t, q), 0.5, 6, {"out": i}, P_init=Pv)
    assert r.finite
    assert r.Q["out"][-1] == pytest.approx(q, rel=1e-4)
    assert r.P["out"][-1] == pytest.approx(Pv + R * q, rel=1e-3)
    assert net.state["Pc"][0] == pytest.approx(Pv - net.P_ref)


def test_matched_resistance_outlet_absorbs_a_pulse():
    net = Network([Vessel("tube", 1.0, 0.01, 0.01, None, c_in=6.0, c_out=6.0)], mu=0.0, dx=0.0025)  # R1 None: matched
    r = net.simulate(_pulse, 0.4, 1, {"mid": net.node("tube", 0.5)}, P_init=0.0)
    p = r.P["mid"]
    assert np.abs(p[r.t > 0.2]).max() < 0.03 * p[r.t < 0.17].max()


def test_windkessel_needs_both_r2_and_c():
    with pytest.raises(AssertionError):
        Network([Vessel("v", 0.1, 3e-3, 3e-3, None, c_in=6.0, c_out=6.0, R2=1e9)])
