#!/usr/bin/env python3
"""Reproduce the illustrative Test 3 information calculations.

Companion to main_seq2.tex, Biological Invariants Under a Singular-Set-Mediated
ECL: A SEQUENTION Account of Canalization, Convergence, and Operational Time.

These are deterministic calculations for declared model assumptions, not
biological observations. Normal-approximation sample sizes are planning values,
not exact minimal sizes or demonstrated experimental feasibility. All parameters
are treated as known; no calibration error, clustering, multiplicity, censoring,
or random exposure correction is included.

Usage (Python 3.10 or later, standard library only):
    python test3_feasibility.py
    python test3_feasibility.py --json
    python test3_feasibility.py --check

Counts m are average expected substitutions PER DIRECTION. Counts R are
independent invasion trials PER ORIENTATION. Sum the two directions/orientations
for a total experimental count. The mutation waiting time is not assigned a
numerical value because no empirical mutation intensity is supplied.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from statistics import NormalDist
from typing import Sequence


@dataclass(frozen=True)
class Scenario:
    n: int
    s: float
    y: float = 0.2

    def validate(self) -> None:
        if isinstance(self.n, bool) or not isinstance(self.n, int) or self.n < 2:
            raise ValueError("n must be an integer >= 2")
        if not math.isfinite(self.s) or self.s == 0:
            raise ValueError("planning scenarios require finite, nonzero s")
        if not math.isfinite(self.y) or self.y <= 0:
            raise ValueError("y = g/a_dagger must be finite and positive")


@dataclass(frozen=True)
class Plan:
    n: int
    s: float
    y: float
    rho_E_plus: float
    rho_E_minus: float
    rho_Q_plus: float
    rho_Q_minus: float
    Delta_E: float
    Delta_Q: float
    delta_E: float
    delta_Q: float
    m_E_vs_Q_per_direction: int
    R_E_vs_Q_per_orientation: int
    R_zero_vs_Q_per_orientation: int


def moran_fixation(n: int, s: float) -> float:
    """Exact single-founder Moran probability for log relative fitness s."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if not math.isfinite(s):
        raise ValueError("s must be finite")
    if s == 0.0:
        return 1.0 / n
    a = abs(s)
    positive = math.expm1(-a) / math.expm1(-n * a)
    return positive if s > 0 else positive * math.exp(-(n - 1) * a)


def mu_zero(y: float) -> float:
    """Selected unregularized mobility y/sqrt(1+y^2)."""
    if not math.isfinite(y) or y < 0:
        raise ValueError("y must be finite and nonnegative")
    return y / math.hypot(1.0, y)


def planning_value(v0: float, v1: float, gap: float,
                   alpha: float, power: float) -> int:
    """Equal-allocation Gaussian planning expression, rounded upward."""
    if not 0 < alpha < 1 or not 0.5 < power < 1:
        raise ValueError("require 0 < alpha < 1 and 0.5 < power < 1")
    if gap == 0 or min(v0, v1) < 0:
        raise ValueError("nonzero gap and nonnegative variances are required")
    normal = NormalDist()
    z0 = normal.inv_cdf(1 - alpha / 2)
    z1 = normal.inv_cdf(power)
    return math.ceil((z0 * math.sqrt(v0) + z1 * math.sqrt(v1)) ** 2 / gap ** 2)


def calculate(scenario: Scenario, alpha: float = 0.05,
              power: float = 0.80) -> Plan:
    scenario.validate()
    n, s, y = scenario.n, scenario.s, scenario.y
    p, r = moran_fixation(n, s), moran_fixation(n, -s)
    mean = (p + r) / 2
    Delta_E = p - r
    Delta_Q = (n - 1) * s * mu_zero(y) / n
    pQ, rQ = mean + Delta_Q / 2, mean - Delta_Q / 2
    if not (0 < pQ < 1 and 0 < rQ < 1):
        raise ValueError("endpoint completion is inadmissible for this scenario")
    delta_E, delta_Q = Delta_E / (2 * mean), Delta_Q / (2 * mean)
    if not abs(delta_Q) < 1:
        raise ValueError("quadratic relative bias lies outside (-1, 1)")
    V_E = p * (1 - p) + r * (1 - r)
    V_Q = pQ * (1 - pQ) + rQ * (1 - rQ)
    V_zero = 2 * mean * (1 - mean)
    # Delta differences concern endpoint probabilities; delta differences concern
    # relative directional biases. They must not be interchanged.
    R_EQ = planning_value(V_E, V_Q, Delta_E - Delta_Q, alpha, power)
    R_0Q = planning_value(V_zero, V_Q, Delta_Q, alpha, power)
    m_EQ = planning_value((1 - delta_E ** 2) / 2,
                          (1 - delta_Q ** 2) / 2,
                          delta_E - delta_Q, alpha, power)
    return Plan(n, s, y, p, r, pQ, rQ, Delta_E, Delta_Q,
                delta_E, delta_Q, m_EQ, R_EQ, R_0Q)


SCENARIOS = (Scenario(10, 0.020), Scenario(20, 0.010),
             Scenario(50, 0.004), Scenario(20, 0.002))


def checks() -> None:
    """Check identities, probability domains, and manuscript table values."""
    for n in (2, 5, 10, 20, 50):
        assert moran_fixation(n, 0) == 1 / n
        for s in (1e-6, 0.002, 0.01, 0.02):
            p, r = moran_fixation(n, s), moran_fixation(n, -s)
            assert 0 < p < 1 and 0 < r < 1
            assert math.isclose(math.log(p / r), (n - 1) * s,
                                rel_tol=1e-9, abs_tol=1e-12)
            assert math.isclose((p-r)/(p+r), math.tanh((n-1)*s/2),
                                rel_tol=1e-8, abs_tol=1e-12)
        h = 1e-6
        derivative = (moran_fixation(n, h)-moran_fixation(n, -h))/(2*h)
        assert math.isclose(derivative, (n-1)/(2*n), rel_tol=1e-7)
    expected = ((750, 6763, 113614), (673, 12821, 215390),
                (633, 31080, 522186), (16823, 319658, 5370865))
    for scenario, counts in zip(SCENARIOS, expected):
        plan = calculate(scenario)
        assert (plan.m_E_vs_Q_per_direction, plan.R_E_vs_Q_per_orientation,
                plan.R_zero_vs_Q_per_orientation) == counts
        assert math.isclose(plan.rho_Q_plus + plan.rho_Q_minus,
                            plan.rho_E_plus + plan.rho_E_minus)
    # Poisson signal-to-noise convention: expected total is 2*m.
    delta = 0.02
    assert math.isclose((1-delta**2)/(2*delta**2), 1249.5)
    assert math.isclose(9*(1-delta**2)/(2*delta**2), 11245.5)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="print full values as JSON")
    parser.add_argument("--check", action="store_true", help="run identity and table checks only")
    args = parser.parse_args(argv)
    checks()
    if args.check:
        print("All identity, probability-domain, and table-reproduction checks passed.")
        return 0
    plans = [calculate(scenario) for scenario in SCENARIOS]
    if args.json:
        print(json.dumps({"kind": "illustrative analytical planning calculations",
                          "alpha_two_sided": 0.05, "nominal_power": 0.80,
                          "counts_are_per_direction_or_orientation": True,
                          "nuisance_uncertainty_included": False,
                          "scenarios": [asdict(plan) for plan in plans]}, indent=2))
    else:
        print("Illustrative planning values (not empirical observations)")
        print(" n      s      y   delta_E   delta_Q       m_EQ        R_EQ         R_0Q")
        for p in plans:
            print(f"{p.n:2d}  {p.s:.4f}  {p.y:.3f}  {p.delta_E:.5f}  {p.delta_Q:.5f}"
                  f"  {p.m_E_vs_Q_per_direction:9,d}"
                  f"  {p.R_E_vs_Q_per_orientation:10,d}"
                  f"  {p.R_zero_vs_Q_per_orientation:11,d}")
        print("m: mean substitution count/direction; R: invasion trials/orientation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
