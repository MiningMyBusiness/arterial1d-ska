"""arterial1d: one-dimensional blood flow in elastic arterial trees.

An independent implementation (two-step MacCormack finite differences) that reads networks in
openBF's YAML format. See README.md for the scheme, validation and limitations.
"""
from importlib.resources import files as _files

from .openbf import UnsupportedFeature, load_inlet, load_network
from .solver import MMHG, Network, Result, Vessel, wall_thickness

__version__ = "0.1.0.dev0"


def adan56_path():
    """Path to the bundled ADAN56 network file (from openBF, Apache-2.0)."""
    return _files(__package__) / "data" / "adan56" / "adan56.yaml"


__all__ = ["MMHG", "Network", "Result", "Vessel", "wall_thickness", "load_network", "load_inlet",
           "UnsupportedFeature", "adan56_path", "__version__"]
