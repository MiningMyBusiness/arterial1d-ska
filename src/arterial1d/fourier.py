"""Inflow waveforms given as formulas, as in the Boileau et al. 2015 supplementary `inflow.txt`
files (cases 2-5): a Fourier series in t with constants PI and T, e.g.

    # period T=1.1
    inflow=(6.5+3.294*sin(2*PI*t/T-0.023974)+...)

The expression is evaluated with a small whitelist interpreter (numbers, t, T, PI, sin, cos,
+ - * / and unary minus), never with eval().
"""
from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Callable, Tuple

import numpy as np

_FUNCS = {"sin": np.sin, "cos": np.cos}
_BIN = {ast.Add: np.add, ast.Sub: np.subtract, ast.Mult: np.multiply, ast.Div: np.divide}


def _compile(node, names):
    if isinstance(node, ast.Expression):
        return _compile(node.body, names)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        v = float(node.value)
        return lambda env: v
    if isinstance(node, ast.Name) and node.id in names:
        k = node.id
        return lambda env: env[k]
    if isinstance(node, ast.BinOp) and type(node.op) in _BIN:
        f, a, b = _BIN[type(node.op)], _compile(node.left, names), _compile(node.right, names)
        return lambda env: f(a(env), b(env))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        a, sgn = _compile(node.operand, names), (-1.0 if isinstance(node.op, ast.USub) else 1.0)
        return lambda env: sgn * a(env)
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in _FUNCS
            and len(node.args) == 1 and not node.keywords):
        f, a = _FUNCS[node.func.id], _compile(node.args[0], names)
        return lambda env: f(a(env))
    raise ValueError(f"unsupported expression element: {ast.dump(node)[:80]}")


def read_inflow_formula(path) -> Tuple[Callable[[np.ndarray], np.ndarray], float, str]:
    """Returns (q(t) in the file's units, period T [s], the file's comment lines)."""
    text = Path(path).read_text()
    comments = "\n".join(l for l in text.splitlines() if l.lstrip().startswith("#"))
    m = re.search(r"period\s*T\s*=\s*([0-9.eE+-]+)", comments)
    if not m:
        raise ValueError(f"{path}: no 'period T=' line")
    period = float(m.group(1))
    m = re.search(r"inflow\s*=\s*(.+)", "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#")), re.S)
    if not m:
        raise ValueError(f"{path}: no 'inflow=' expression")
    expr = re.sub(r"\s+", "", m.group(1))
    f = _compile(ast.parse(expr, mode="eval"), {"t", "T", "PI"})
    return (lambda t: np.asarray(f({"t": np.asarray(t, float), "T": period, "PI": np.pi}), float)
            * np.ones_like(np.asarray(t, float))), period, comments
