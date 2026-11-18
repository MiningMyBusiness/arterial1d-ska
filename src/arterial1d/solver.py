"""One-dimensional blood flow in a tree of elastic arteries.

Equations, per vessel (Sherwin et al. 2003; Alastruey et al. 2012):
    dA/dt + dQ/dx = 0
    dQ/dt + d(alpha Q^2/A)/dx + (A/rho) dP/dx = -K_R Q/A,      K_R = 2 (gamma + 2) pi mu / rho
    P = P_ext + beta(x) (sqrt(A) - sqrt(A0(x))),               c = sqrt(beta/(2 rho)) A^(1/4)

gamma is the velocity-profile exponent: gamma = 2 (parabolic, K_R = 8 pi mu / rho) is the
default; gamma = 9 gives the flatter profile often used for large arteries (K_R = 22 pi mu / rho).

Numerics: two-step MacCormack (second order), with the pressure term kept as (A/rho) dP/dx so
that a uniform-pressure rest state is preserved exactly on tapered vessels. Boundaries use the
characteristic compatibility relations in pressure-velocity form, dP +/- rho c du, each reduced
to one scalar Newton solve for the boundary pressure:
    inlet     prescribed flow Q_in(t)
    junction  mass conservation + continuity of total pressure P + rho u^2/2 (static pressure
              optional), one parent and up to 3 daughters
    outlet    three-element Windkessel:  Q = (P - P_c)/R1,  C dP_c/dt = Q - (P_c - P_ven)/R2,
              or a pure resistance Q = (P - P_ven)/R1 (no C, no R2); with R1 matched to the
              characteristic impedance this is a (linearly) non-reflecting outlet
The time step is fixed for a run (see Network.dt_max).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import warnings

import numpy as np
from numba import njit

MMHG = 133.322
MAX_DAUGHTERS = 3


def wall_thickness(r):
    """Empirical wall thickness [m] for radius r [m] (the law used with ADAN56 in Boileau et al. 2015 and openBF)."""
    return r * (0.2802 * np.exp(-505.3 * r) + 0.1324 * np.exp(-11.14 * r))


@dataclass
class Vessel:
    name: str
    length: float                 # m
    r_in: float                   # unstressed radius at the upstream end, m
    r_out: float                  # at the downstream end, m
    parent: Optional[str] = None  # None for the root vessel
    # stiffness: either reference wave speeds at both ends, or a Young's modulus
    c_in: Optional[float] = None  # wave speed at A0, upstream end, m/s
    c_out: Optional[float] = None
    E: Optional[float] = None     # Pa; c0 = sqrt(2 E h / (3 rho r))
    h: Optional[float] = None     # wall thickness, m (with E); None = the empirical thickness law
    taper: str = "exp"            # radius interpolation along the vessel: "exp" or "linear"
    # terminal condition (only for vessels with no daughters): three-element Windkessel R1-C-R2,
    # or a pure resistance R1 when R2 and C are both None
    R1: Optional[float] = None    # Pa s/m^3; None = matched to the characteristic impedance
    R2: Optional[float] = None
    C: Optional[float] = None     # m^3/Pa


@dataclass
class Result:
    t: np.ndarray                         # s, recorded samples
    P: Dict[str, np.ndarray]              # Pa, per probe
    A: Dict[str, np.ndarray]              # m^2
    Q: Dict[str, np.ndarray]              # m^3/s
    dt: float
    volume_in: np.ndarray                 # m^3 per cycle
    volume_out: np.ndarray                # m^3 per cycle, summed over outlets
    volume_out_by_outlet: Dict[str, np.ndarray]
    stored_volume: np.ndarray             # m^3 in the 1D vessels at the end of each cycle
    finite: bool
    max_courant: float = float("nan")     # max over nodes and steps of (|u| + c) dt / dx (stable below 1)


@njit(cache=True)
def _area(P, A0, beta):
    s = np.sqrt(A0) + P / beta
    return s * s


@njit(cache=True)
def _run(A, Q, A0, beta, start, end, dxv, parent, daughters, n_daughters, is_terminal,
         R1, R2, C, Pc, P_ven, q_in, dt, rho, K_R, alpha, probes, rec_every, steps_per_cycle,
         total_pressure, cb, inv_dxn):
    n_v = len(start)
    n_steps = len(q_in)
    n_rec = n_steps // rec_every
    n_cyc = n_steps // steps_per_cycle
    recP = np.zeros((len(probes), n_rec))
    recA = np.zeros((len(probes), n_rec))
    recQ = np.zeros((len(probes), n_rec))
    vol_in = np.zeros(n_cyc)
    vol_out = np.zeros((n_v, n_cyc))
    stored = np.zeros(n_cyc)
    As = np.empty_like(A)
    Qs = np.empty_like(Q)
    An = np.empty_like(A)
    Qn = np.empty_like(Q)
    sqA0 = np.sqrt(A0)
    off_d = np.zeros(MAX_DAUGHTERS)
    ok = True
    courant = 0.0
    for k in range(n_steps):
        P = beta * (np.sqrt(A) - sqA0)
        for i in range(len(A)):                 # (|u| + c) dt / dx, c = sqrt(beta/(2 rho)) A^(1/4)
            cr = (abs(Q[i] / A[i]) + cb[i] * np.sqrt(np.sqrt(A[i]))) * dt * inv_dxn[i]
            if cr > courant:
                courant = cr
        # ---- interior: MacCormack, direction alternating each step -------------------------
        fwd = (k % 2 == 0)
        for v in range(n_v):
            s, e = start[v], end[v]
            r = dt / dxv[v]
            for i in range(s, e):
                j = i + 1 if fwd else i - 1
                if j < s or j >= e:
                    As[i] = A[i]
                    Qs[i] = Q[i]
                    continue
                sg = 1.0 if fwd else -1.0
                dF = (alpha * Q[j] * Q[j] / A[j] - alpha * Q[i] * Q[i] / A[i]) * sg
                As[i] = A[i] - r * (Q[j] - Q[i]) * sg
                Qs[i] = Q[i] - r * (dF + A[i] / rho * (P[j] - P[i]) * sg) - dt * K_R * Q[i] / A[i]
            for i in range(s + 1, e - 1):
                j = i - 1 if fwd else i + 1
                sg = 1.0 if fwd else -1.0
                Pi = beta[i] * (np.sqrt(As[i]) - sqA0[i])
                Pj = beta[j] * (np.sqrt(As[j]) - sqA0[j])
                dF = (alpha * Qs[i] * Qs[i] / As[i] - alpha * Qs[j] * Qs[j] / As[j]) * sg
                An[i] = 0.5 * (A[i] + As[i] - r * (Qs[i] - Qs[j]) * sg)
                Qn[i] = 0.5 * (Q[i] + Qs[i] - r * (dF + As[i] / rho * (Pi - Pj) * sg)
                               - dt * K_R * Qs[i] / As[i])
        # ---- boundaries: compatibility relations from time level n -------------------------
        # R_up[v]: backward-running relation at the upstream end,  P - rho c u = R_up
        # R_dn[v]: forward-running relation at the downstream end, P + rho c u = R_dn
        for v in range(n_v):
            s, e = start[v], end[v]
            # upstream end of v
            c = np.sqrt(beta[s] / (2 * rho)) * A[s] ** 0.25
            u = Q[s] / A[s]
            xi = (c - u) * dt / dxv[v]
            Pf = P[s] + xi * (P[s + 1] - P[s])
            uf = u + xi * (Q[s + 1] / A[s + 1] - u)
            Af = A[s] + xi * (A[s + 1] - A[s])
            # taper term: along both characteristic families d(P -/+ rho c u)/dt carries
            # G = -rho c^2 u (dA/dx at constant P) / A, zero at rest and in a uniform vessel
            G = -rho * c * c * uf * (_area(Pf, A0[s + 1], beta[s + 1]) - _area(Pf, A0[s], beta[s])) / (dxv[v] * Af)
            Rup = Pf - rho * c * uf + dt * rho * c * K_R * uf / Af + dt * G
            zc_up = rho * c
            # downstream end of v
            m = e - 1
            c2 = np.sqrt(beta[m] / (2 * rho)) * A[m] ** 0.25
            u2 = Q[m] / A[m]
            xi2 = (c2 + u2) * dt / dxv[v]
            Pf2 = P[m] - xi2 * (P[m] - P[m - 1])
            uf2 = u2 - xi2 * (u2 - Q[m - 1] / A[m - 1])
            Af2 = A[m] - xi2 * (A[m] - A[m - 1])
            G2 = -rho * c2 * c2 * uf2 * (_area(Pf2, A0[m], beta[m]) - _area(Pf2, A0[m - 1], beta[m - 1])) / (dxv[v] * Af2)
            Rdn = Pf2 + rho * c2 * uf2 - dt * rho * c2 * K_R * uf2 / Af2 + dt * G2
            zc_dn = rho * c2
            # stash in the predictor arrays' boundary slots (not otherwise used past this point)
            As[s] = Rup
            Qs[s] = zc_up
            As[m] = Rdn
            Qs[m] = zc_dn
        for v in range(n_v):
            s, e = start[v], end[v]
            m = e - 1
            if parent[v] < 0:
                # inlet: A(P) (P - Rup)/(rho c) = Q_in
                Rup, zc = As[s], Qs[s]
                Pb = P[s]
                for _ in range(8):
                    Ab = _area(Pb, A0[s], beta[s])
                    f = Ab * (Pb - Rup) / zc - q_in[k]
                    df = (2 * np.sqrt(Ab) / beta[s] * (Pb - Rup) + Ab) / zc
                    step = f / df
                    Pb -= step
                    if abs(step) < 1e-9 * (abs(Pb) + 1.0):
                        break
                An[s] = _area(Pb, A0[s], beta[s])
                Qn[s] = q_in[k]
            if is_terminal[v]:
                # outlet: A(P) (Rdn - P)/(rho c) = (P - Pc)/R1
                Rdn, zc = As[m], Qs[m]
                Pb = P[m]
                for _ in range(8):
                    Ab = _area(Pb, A0[m], beta[m])
                    f = Ab * (Rdn - Pb) / zc - (Pb - Pc[v]) / R1[v]
                    df = (2 * np.sqrt(Ab) / beta[m] * (Rdn - Pb) - Ab) / zc - 1.0 / R1[v]
                    step = f / df
                    Pb -= step
                    if abs(step) < 1e-9 * (abs(Pb) + 1.0):
                        break
                An[m] = _area(Pb, A0[m], beta[m])
                Qn[m] = (Pb - Pc[v]) / R1[v]
                if C[v] > 0.0:                      # Windkessel; otherwise Pc stays at P_ven
                    Pc[v] += dt * (Qn[m] - (Pc[v] - P_ven) / R2[v]) / C[v]
            else:
                # junction at the downstream end of v: mass conserved and, with total_pressure,
                # P + rho u^2 / 2 continuous (static pressure otherwise). The dynamic-pressure
                # offsets are lagged and refined in a short outer iteration.
                Rdn, zc = As[m], Qs[m]
                Pb = P[m]
                off_p = 0.0
                nd = n_daughters[v]
                for q in range(nd):
                    off_d[q] = 0.0
                if total_pressure:
                    up = Q[m] / A[m]
                    off_p = -0.5 * rho * up * up
                    Pb = P[m] - off_p
                    for q in range(nd):
                        sd = start[daughters[v, q]]
                        ud = Q[sd] / A[sd]
                        off_d[q] = -0.5 * rho * ud * ud
                for outer in range(4 if total_pressure else 1):
                    for _ in range(10):
                        Pp = Pb + off_p
                        Ab = _area(Pp, A0[m], beta[m])
                        f = Ab * (Rdn - Pp) / zc
                        df = (2 * np.sqrt(Ab) / beta[m] * (Rdn - Pp) - Ab) / zc
                        for q in range(nd):
                            d = daughters[v, q]
                            sd = start[d]
                            Pd = Pb + off_d[q]
                            Ad = _area(Pd, A0[sd], beta[sd])
                            f -= Ad * (Pd - As[sd]) / Qs[sd]
                            df -= (2 * np.sqrt(Ad) / beta[sd] * (Pd - As[sd]) + Ad) / Qs[sd]
                        step = f / df
                        Pb -= step
                        if abs(step) < 1e-9 * (abs(Pb) + 1.0):
                            break
                    if total_pressure:
                        up = (Rdn - (Pb + off_p)) / zc
                        off_p = -0.5 * rho * up * up
                        for q in range(nd):
                            sd = start[daughters[v, q]]
                            ud = (Pb + off_d[q] - As[sd]) / Qs[sd]
                            off_d[q] = -0.5 * rho * ud * ud
                # final values: daughters from their own relations, parent flow from mass conservation
                qsum = 0.0
                for q in range(nd):
                    d = daughters[v, q]
                    sd = start[d]
                    Pd = Pb + off_d[q]
                    Ad = _area(Pd, A0[sd], beta[sd])
                    An[sd] = Ad
                    Qn[sd] = Ad * (Pd - As[sd]) / Qs[sd]
                    qsum += Qn[sd]
                An[m] = _area(Pb + off_p, A0[m], beta[m])
                Qn[m] = qsum
        # ---- commit, bookkeeping -----------------------------------------------------------
        for i in range(len(A)):
            A[i] = An[i]
            Q[i] = Qn[i]
        cyc = k // steps_per_cycle
        if cyc < n_cyc:
            vol_in[cyc] += q_in[k] * dt
            for v in range(n_v):
                if is_terminal[v]:
                    vol_out[v, cyc] += Q[end[v] - 1] * dt
            if (k + 1) % steps_per_cycle == 0:
                tot = 0.0
                for v in range(n_v):
                    s, e = start[v], end[v]
                    tot += dxv[v] * (np.sum(A[s:e]) - 0.5 * (A[s] + A[e - 1]))
                stored[cyc] = tot
        if (k + 1) % rec_every == 0:
            j = (k + 1) // rec_every - 1
            if j < n_rec:
                for p in range(len(probes)):
                    i = probes[p]
                    recP[p, j] = beta[i] * (np.sqrt(A[i]) - sqA0[i])
                    recA[p, j] = A[i]
                    recQ[p, j] = Q[i]
        if not (np.isfinite(A[0]) and np.isfinite(Q[0]) and np.isfinite(A[len(A) - 1])) or np.min(A) <= 0.0:
            ok = False
            break
    return recP, recA, recQ, vol_in, vol_out, stored, ok, courant


class Network:
    """A tree of tapered vessels with a flow inlet at the root and Windkessel outlets."""

    def __init__(self, vessels: List[Vessel], rho=1060.0, mu=0.004, alpha=1.0, friction_coeff_pi=8.0,
                 dx=0.0025, cfl=0.9, P_ven=0.0, P_ref=0.0, min_nodes=6, c_headroom=1.25,
                 total_pressure=True):
        """vessels: root first, every vessel after its parent. Units SI.

        friction_coeff_pi: K_R in units of pi mu / rho, 2 (gamma + 2) for profile exponent gamma.
        dx: target node spacing [m]; each vessel gets at least min_nodes nodes spanning its length.
            The single time step is set by the smallest spacing in the network (often a short
            vessel), so the Courant number elsewhere is lower; see README (case 1) for its effect.
        cfl, c_headroom: the fixed time step is cfl * dx / (c_headroom * c0), with c0 the wave
            speed at the reference area, so that the CFL condition holds while c rises with
            pressure (c ~ A^(1/4)) and flow speed adds to it. The defaults give a Courant number of
            about 0.72 at rest; on ADAN56 the maximum reached is 1.28x the resting value (0.92).
            Lower Courant numbers damp short waveform features more (Boileau case 1; README).
            Result.max_courant reports the maximum reached; above 1 simulate() warns.
        total_pressure: junctions conserve P + rho u^2 / 2 (True) or static pressure P (False).
        P_ref is the pressure at which the vessels have their reference areas A0 (0 for
        unstressed radii, the diastolic pressure for ADAN56). Pressures passed in and reported
        are absolute; the solver works with P - P_ref."""
        self.vessels = vessels
        self.total_pressure = total_pressure
        self.rho, self.mu, self.alpha, self.P_ven, self.P_ref = rho, mu, alpha, P_ven, P_ref
        self.K_R = friction_coeff_pi * np.pi * mu / rho
        names = [v.name for v in vessels]
        assert len(set(names)) == len(names)
        self.index = {n: i for i, n in enumerate(names)}
        roots = [v for v in vessels if v.parent is None]
        assert len(roots) == 1 and vessels[0].parent is None, "first vessel must be the single root"
        n_v = len(vessels)
        self.start = np.zeros(n_v, np.int64)
        self.end = np.zeros(n_v, np.int64)
        self.dxv = np.zeros(n_v)
        A0, beta, xs = [], [], []
        pos = 0
        for i, v in enumerate(vessels):
            nx = max(int(round(v.length / dx)) + 1, min_nodes)
            x = np.linspace(0.0, v.length, nx)     # the grid always spans the vessel length
            if v.taper == "linear":
                r0 = v.r_in + (v.r_out - v.r_in) * x / v.length
            else:
                r0 = v.r_in * (v.r_out / v.r_in) ** (x / v.length)
            if v.E is not None:
                h = wall_thickness(r0) if v.h is None else v.h
                c0 = np.sqrt(2.0 * v.E * h / (3.0 * rho * r0))
            else:
                c0 = v.c_in * (v.c_out / v.c_in) ** (x / v.length)
            a0 = np.pi * r0 ** 2
            A0.append(a0)
            beta.append(2.0 * rho * c0 ** 2 / np.sqrt(a0))   # from c0 = sqrt(beta/(2 rho)) A0^(1/4)
            xs.append(x)
            self.start[i], self.end[i], self.dxv[i] = pos, pos + nx, x[1] - x[0]
            pos += nx
        self.A0, self.beta, self.x = np.concatenate(A0), np.concatenate(beta), xs
        self.parent = np.array([-1 if v.parent is None else self.index[v.parent] for v in vessels], np.int64)
        self.daughters = -np.ones((n_v, MAX_DAUGHTERS), np.int64)
        self.n_daughters = np.zeros(n_v, np.int64)
        for i, p in enumerate(self.parent):
            if p >= 0:
                assert self.n_daughters[p] < MAX_DAUGHTERS, f"{vessels[p].name}: too many daughters"
                self.daughters[p, self.n_daughters[p]] = i
                self.n_daughters[p] += 1
        self.is_terminal = self.n_daughters == 0
        self.R1, self.R2, self.C = np.ones(n_v), np.ones(n_v), np.ones(n_v)
        for i, v in enumerate(vessels):
            if self.is_terminal[i]:
                assert (v.R2 is None) == (v.C is None), f"{v.name}: give both R2 and C, or neither (pure resistance)"
                m = self.end[i] - 1
                z0 = np.sqrt(self.beta[m] / (2 * rho)) * self.A0[m] ** 0.25 * rho / self.A0[m]
                self.R1[i] = z0 if v.R1 is None else v.R1
                if v.C is None:
                    self.R2[i], self.C[i] = np.inf, 0.0      # pure resistance
                else:
                    assert v.C > 0 and v.R2 > 0, f"{v.name}: R2 and C must be positive"
                    self.R2[i], self.C[i] = v.R2, v.C
        # fixed time step from the stiffest state expected: c grows as A^(1/4), allow headroom
        c0 = np.sqrt(self.beta / (2 * rho)) * self.A0 ** 0.25
        self.dt_max = cfl * np.min(np.repeat(self.dxv, self.end - self.start) / (c_headroom * c0))

    def node(self, vessel: str, frac: float = 1.0) -> int:
        """Flat node index at a fractional position along a vessel (0 = upstream, 1 = downstream)."""
        i = self.index[vessel]
        n = self.end[i] - self.start[i]
        return int(self.start[i] + min(n - 1, max(0, round(frac * (n - 1)))))

    def wave_speed0(self) -> np.ndarray:
        return np.sqrt(self.beta / (2 * self.rho)) * self.A0 ** 0.25

    def simulate(self, q_in: Callable[[np.ndarray], np.ndarray], period: float, n_cycles: int,
                 probes: Dict[str, int], P_init: float = 80 * MMHG, record_hz: float = 1000.0,
                 state: Optional[dict] = None) -> Result:
        """Run n_cycles of period `period` [s]. q_in maps time [s] to inlet flow [m^3/s].

        The step is the largest dt <= dt_max that divides both the period and the recording
        interval, so every cycle has the same number of steps and samples fall on step ends.
        """
        rec_dt = 1.0 / record_hz
        n_rec_per_cycle = int(round(period / rec_dt))
        assert abs(n_rec_per_cycle * rec_dt - period) < 1e-9, "period must be a multiple of 1/record_hz"
        rec_every = int(np.ceil(rec_dt / self.dt_max))
        dt = rec_dt / rec_every
        spc = n_rec_per_cycle * rec_every
        tq = (np.arange(n_cycles * spc) + 1) * dt          # flow imposed at the new time level
        qin = np.ascontiguousarray(q_in(tq), dtype=np.float64)
        if state is None:
            A = _area(np.full_like(self.A0, P_init - self.P_ref), self.A0, self.beta)
            Q = np.zeros_like(A)
            Pc = np.full(len(self.vessels), P_init - self.P_ref)
            Pc[self.C == 0.0] = self.P_ven - self.P_ref          # pure-resistance outlets
        else:
            A, Q, Pc = state["A"].copy(), state["Q"].copy(), state["Pc"].copy()
        names = list(probes)
        pidx = np.array([probes[n] for n in names], np.int64)
        cb = np.sqrt(self.beta / (2 * self.rho))
        inv_dxn = np.repeat(1.0 / self.dxv, self.end - self.start)
        recP, recA, recQ, vin, vout, stored, ok, courant = _run(
            A, Q, self.A0, self.beta, self.start, self.end, self.dxv, self.parent, self.daughters,
            self.n_daughters, self.is_terminal, self.R1, self.R2, self.C, Pc, self.P_ven - self.P_ref, qin, dt,
            self.rho, self.K_R, self.alpha, pidx, rec_every, spc, self.total_pressure, cb, inv_dxn)
        if courant > 1.0:
            warnings.warn(f"Courant number reached {courant:.2f} (> 1): the run is unreliable; "
                          f"increase c_headroom or reduce cfl", RuntimeWarning, stacklevel=2)
        self.state = {"A": A, "Q": Q, "Pc": Pc}
        t = (np.arange(recP.shape[1]) + 1) * rec_dt
        return Result(
            t=t, P={n: recP[i] + self.P_ref for i, n in enumerate(names)}, A={n: recA[i] for i, n in enumerate(names)},
            Q={n: recQ[i] for i, n in enumerate(names)}, dt=dt, volume_in=vin,
            volume_out=vout[self.is_terminal].sum(axis=0),
            volume_out_by_outlet={v.name: vout[i] for i, v in enumerate(self.vessels) if self.is_terminal[i]},
            stored_volume=stored, finite=bool(ok), max_courant=float(courant))
