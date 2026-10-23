"""physverify spec for the driven LCR circuit in v2_modern_python/.

    physverify run tests/verify_lcr.py

Equation, from Kirchhoff's voltage law as written in the README:

    L q'' + R q' + q / C = V0 sin(omega t)

with L = 1 H, C = 0.0025 F, omega = 2.5 rad/s, V0 = 0.2 V, as a first-order
system for y = [q, q']. (The code names its variables current and di/dt and its
forcing E_dot; with those names it is the same equation for i, driven by
dE/dt = V0 sin(omega t). The numbers are identical either way.)

The equation is linear with constant coefficients, so it has an exact solution:
a decaying homogeneous part plus the steady state A sin(omega t) + B cos(omega t).

What is checked:
- LCRSimulator.solve_rk4: manufactured solution, reference solution, whole
  trajectory and exact patch test, at the RK4 order 4, for R = 10 ohm.
- The time grid of solve_rk4 ends at t_end after exactly n steps.
- All three damping cases of the README (R = 10, 40, 100 ohm: under-, critically
  and over-damped) in one vectorized call: the largest error against the exact
  solution over 0 to 2 s falls at order 4.
- The README validation (R = 10 ohm, h = 0.01 s, 0 to 5 s): solve_rk4 against the
  exact solution, and scipy_validation.py's own RK4 against SciPy RK45 (the
  README's 9.97e-8), each with a maximum error below 1e-7.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "v2_modern_python"))

from lcr_engine import LCRSimulator  # noqa: E402
from physverify import Check, ManufacturedODE, Result, check_order, ode_checks  # noqa: E402
from physverify.studies import time_grid_check  # noqa: E402

L, C, OMEGA, V0 = 1.0, 0.0025, 2.5, 0.2
R = 10.0
Y0 = np.array([0.0, 0.1])
T_END = 1.0
NS = (20, 40, 80, 160)


def reference_rhs(t, y, r=R):
    q, i = y
    return np.array([i, (V0 * np.sin(OMEGA * t) - r * i - q / C) / L])


def exact(t, r, y0=Y0):
    """Exact q(t) and q'(t) for L q'' + r q' + q/C = V0 sin(omega t)."""
    t = np.asarray(t, dtype=float)
    gamma, w0sq = r / (2 * L), 1.0 / (L * C)
    f = V0 / L
    den = (w0sq - OMEGA**2) ** 2 + (2 * gamma * OMEGA) ** 2
    a, b = f * (w0sq - OMEGA**2) / den, -f * 2 * gamma * OMEGA / den
    qp = a * np.sin(OMEGA * t) + b * np.cos(OMEGA * t)
    ip = OMEGA * (a * np.cos(OMEGA * t) - b * np.sin(OMEGA * t))
    q0, i0 = y0[0] - b, y0[1] - a * OMEGA  # homogeneous part at t = 0
    disc = gamma**2 - w0sq
    e = np.exp(-gamma * t)
    if abs(disc) < 1e-12:  # critical damping: (c1 + c2 t) e^{-gamma t}
        c1, c2 = q0, i0 + gamma * q0
        qh = (c1 + c2 * t) * e
        ih = (c2 - gamma * (c1 + c2 * t)) * e
    elif disc < 0:  # under-damped
        wd = np.sqrt(-disc)
        c1, c2 = q0, (i0 + gamma * q0) / wd
        qh = e * (c1 * np.cos(wd * t) + c2 * np.sin(wd * t))
        ih = -gamma * qh + e * wd * (-c1 * np.sin(wd * t) + c2 * np.cos(wd * t))
    else:  # over-damped: two real decay rates
        s1, s2 = -gamma + np.sqrt(disc), -gamma - np.sqrt(disc)
        c1 = (i0 - s2 * q0) / (s1 - s2)
        c2 = q0 - c1
        qh = c1 * np.exp(s1 * t) + c2 * np.exp(s2 * t)
        ih = c1 * s1 * np.exp(s1 * t) + c2 * s2 * np.exp(s2 * t)
    return np.array([qp + qh, ip + ih])


class _Forced(LCRSimulator):
    """LCRSimulator with an extra term source(t) on the right-hand side."""

    def __init__(self, source):
        super().__init__(L=L, C=C, omega=OMEGA)
        self.source = source

    def vectorized_derivatives(self, t, state, R_array):  # noqa: N803 (name as in the code)
        return super().vectorized_derivatives(t, state, R_array) + np.broadcast_to(self.source(t), (2,)).reshape(2, 1)


def _run(source, y0, t_end, n):
    sim = _Forced(source)
    times, states = sim.solve_rk4((0.0, t_end), t_end / n, np.asarray(y0, dtype=float).reshape(2, 1), np.array([R]))
    return times, states[:, :, 0]


def integrate(source, y0, t_end, n):
    return _run(source, y0, t_end, n)[1][n]


def integrate_trajectory(source, y0, t_end, n):
    return _run(source, y0, t_end, n)[1][: n + 1]


def manufactured(t):
    return np.array([0.01 * np.sin(3.0 * t) + 0.002, 0.03 * np.cos(3.0 * t)])


def damping_cases():
    r_values = np.array([10.0, 40.0, 100.0])
    sim = LCRSimulator(L=L, C=C, omega=OMEGA)
    hs, errs = [0.008, 0.004, 0.002, 0.001], []
    for h in hs:
        t, states = sim.solve_rk4((0.0, 2.0), h, np.tile(Y0.reshape(2, 1), (1, 3)), r_values)
        errs.append(max(float(np.max(np.abs(states[:, :, k] - exact(t, r).T))) for k, r in enumerate(r_values)))
    return check_order(hs, errs, 4.0, scale=float(np.max(np.abs(Y0))))


def readme_validation():
    sim = LCRSimulator(L=L, C=C, omega=OMEGA)
    t, states = sim.solve_rk4((0.0, 5.0), 0.01, Y0.reshape(2, 1), np.array([R]))
    err_exact = float(np.max(np.abs(states[:, 0, 0] - exact(t, R)[0])))

    import scipy_validation  # pylint: disable=import-outside-toplevel
    from scipy.integrate import solve_ivp  # pylint: disable=import-outside-toplevel

    t2, data = scipy_validation.custom_rk4((0, 5.0), 0.01, Y0)
    sol = solve_ivp(scipy_validation.circuit_derivatives, (0, 5.0), Y0, method="RK45", t_eval=t2, rtol=1e-9, atol=1e-9)
    err_scipy = float(np.max(np.abs(data[:, 0] - sol.y[0])))
    err_scipy_exact = float(np.max(np.abs(sol.y[0] - exact(t2, R)[0])))
    text = (f"solve_rk4 vs exact {err_exact:.1e}; scipy_validation RK4 vs SciPy {err_scipy:.2e} "
            f"(SciPy's own error vs exact {err_scipy_exact:.1e})")
    if max(err_exact, err_scipy) < 1e-7:
        return Result("PASS", text)
    return Result("FAIL", text + ", above 1e-7")


def checks():
    out = ode_checks(
        "LCRSimulator.solve_rk4",
        integrate,
        rhs_ref=reference_rhs,
        mms=ManufacturedODE(manufactured, reference_rhs),
        y0_free=Y0,
        t_end=T_END,
        ns=NS,
        order=4.0,
        patch=(np.array([0.01, -0.02]), np.array([0.05, 0.03])),
        integrate_trajectory=integrate_trajectory,
    )
    return out + [
        Check("limit", "solve_rk4 time grid ends at t_end",
              lambda: time_grid_check(lambda n: _run(lambda t: np.zeros(2), Y0, T_END, n)[0], 0.0, T_END,
                                      (3, 4, 7, 8, 10, 16, 20, 100))),
        Check("reference", "solve_rk4 vs exact solution, R = 10, 40, 100 ohm (order 4)", damping_cases),
        Check("reference", "README validation, R = 10 ohm, h = 0.01 s", readme_validation),
    ]
