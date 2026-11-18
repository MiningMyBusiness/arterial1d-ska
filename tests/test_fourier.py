import numpy as np
import pytest

from arterial1d.fourier import read_inflow_formula


def test_reads_a_fourier_series(tmp_path):
    p = tmp_path / "inflow.txt"
    p.write_text("# PI=3.14\n# period T=0.8\n# time t\n\ninflow=1e6*(2.0+0.5*sin(2*PI*t/T-0.1)-0.25*cos(4*PI*t/T+0.2))/1000/60\n")
    f, T, _ = read_inflow_formula(p)
    t = np.linspace(0, 1, 7)
    want = 1e6 * (2.0 + 0.5 * np.sin(2 * np.pi * t / 0.8 - 0.1) - 0.25 * np.cos(4 * np.pi * t / 0.8 + 0.2)) / 1000 / 60
    assert T == 0.8 and np.allclose(f(t), want, rtol=1e-14)


def test_refuses_anything_but_arithmetic(tmp_path):
    p = tmp_path / "inflow.txt"
    p.write_text("# period T=1\ninflow=__import__('os').getcwd()\n")
    with pytest.raises(ValueError):
        read_inflow_formula(p)
