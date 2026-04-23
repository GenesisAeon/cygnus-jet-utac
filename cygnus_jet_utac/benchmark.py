"""Validation against Prabu et al. (2026, Nature Astronomy) measurements.

Runs the CygnusJetUTAC system and scores each observable against the
published reference values. Returns a BenchmarkReport with pass/fail,
% deviation, and an overall score (0–100).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from cygnus_jet_utac.system import CygnusJetUTAC

# (target_value, tolerance_fraction)
PRABU_2026_TARGETS: dict[str, tuple[float, float]] = {
    "jet_power_W":           (3.846e30, 0.15),   # 10 000 L☉ = 3.846e30 W
    "jet_velocity_c":        (0.50,     0.10),
    "accretion_efficiency":  (0.10,     0.05),
    "jet_extent_ly":         (16.0,     0.20),
    "orbital_period_days":   (5.6,      0.01),
    "dance_events_per_year": (2.0,      0.50),
}


@dataclass
class ObservableResult:
    """Single-observable benchmark result.

    Args:
        name: Observable name.
        measured: Value produced by the simulation.
        target: Reference value from Prabu et al. 2026.
        tolerance: Acceptable fractional deviation.
        passed: Whether |deviation| ≤ tolerance.
        deviation_frac: Signed fractional deviation (measured−target)/target.
    """
    name: str
    measured: float
    target: float
    tolerance: float
    passed: bool
    deviation_frac: float
    deviation_pct: float


@dataclass
class BenchmarkReport:
    """Full benchmark report comparing simulation to Prabu et al. 2026.

    Args:
        results: Per-observable results.
        score: Overall score 0–100.
        n_passed: Number of observables that passed.
        n_total: Total number of observables tested.
        gamma_jet: Calibrated Γ_jet value.
        sigma_phi_satisfied: Whether Frame Principle holds.
    """
    results: list[ObservableResult]
    score: float
    n_passed: int
    n_total: int
    gamma_jet: float
    sigma_phi_satisfied: bool
    latex_table: str = field(default="", repr=False)
    zenodo_json: dict = field(default_factory=dict, repr=False)

    def __str__(self) -> str:
        lines = [
            "=" * 64,
            "  BENCHMARK REPORT — Prabu et al. (2026, Nature Astronomy)",
            "=" * 64,
            f"  Overall score: {self.score:.1f}/100  "
            f"({self.n_passed}/{self.n_total} observables passed)",
            f"  Γ_jet = {self.gamma_jet:.4f}  "
            f"| Frame Principle: {'✓' if self.sigma_phi_satisfied else '✗'}",
            "",
            f"  {'Observable':<26} {'Sim':>12} {'Target':>12} {'Dev %':>7} {'Pass':>5}",
            "  " + "-" * 60,
        ]
        for r in self.results:
            tick = "✓" if r.passed else "✗"
            lines.append(
                f"  {r.name:<26} {r.measured:>12.4g} {r.target:>12.4g} "
                f"{r.deviation_pct:>+7.1f}%  {tick}"
            )
        lines.append("=" * 64)
        return "\n".join(lines)


def run_benchmark(system: "CygnusJetUTAC") -> dict:
    """Run full validation of a CygnusJetUTAC instance against Prabu 2026.

    If the system has no results (run_cycle not yet called), runs an 18-year
    cycle automatically.

    Args:
        system: CygnusJetUTAC instance (already run or not).

    Returns:
        BenchmarkReport serialized as a dictionary.
    """
    if not system._results:
        system.run_cycle(duration_years=18.0)

    res = system._results

    # ── Map simulation outputs to benchmark keys ──────────────────────────────
    measured: dict[str, float] = {
        "jet_power_W":           res.get("jet_power_W", 0.0),
        "jet_velocity_c":        res.get("jet_velocity_c", 0.0),
        "accretion_efficiency":  res.get("accretion_efficiency", 0.0),
        "jet_extent_ly":         res.get("jet_extent_ly", 0.0),
        "orbital_period_days":   res.get("orbital_period_days", 0.0),
        "dance_events_per_year": res.get("dance_events_per_year", 0.0),
    }

    obs_results: list[ObservableResult] = []
    n_passed = 0

    for name, (target, tol) in PRABU_2026_TARGETS.items():
        sim_val = measured.get(name, 0.0)
        if target != 0.0:
            dev_frac = (sim_val - target) / target
        else:
            dev_frac = float("inf")
        dev_pct = dev_frac * 100.0
        passed = abs(dev_frac) <= tol
        if passed:
            n_passed += 1
        obs_results.append(ObservableResult(
            name=name,
            measured=sim_val,
            target=target,
            tolerance=tol,
            passed=passed,
            deviation_frac=dev_frac,
            deviation_pct=dev_pct,
        ))

    n_total = len(obs_results)
    # Score: base 100 for all pass; deduct for failures weighted by relative deviation
    score = 100.0 * n_passed / n_total
    # Partial credit: for failed items, give partial credit if within 2× tolerance
    for r in obs_results:
        if not r.passed:
            over = abs(r.deviation_frac) / r.tolerance  # how many tolerances over
            partial = max(0.0, 1.0 - (over - 1.0))  # 0 at 2× tolerance
            score += (100.0 / n_total) * partial * 0.5

    score = min(100.0, score)

    # ── LaTeX table ───────────────────────────────────────────────────────────
    latex_lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Benchmark against Prabu et al. (2026, Nature Astronomy)}",
        r"\begin{tabular}{lrrrc}",
        r"\hline",
        r"Observable & Simulation & Target & Dev.\,(\%) & Pass \\",
        r"\hline",
    ]
    for r in obs_results:
        tick = r"\checkmark" if r.passed else r"\times"
        safe_name = r.name.replace("_", r"\_")
        latex_lines.append(
            f"{safe_name} & {r.measured:.3g} & "
            f"{r.target:.3g} & {r.deviation_pct:+.1f} & ${tick}$ \\\\"
        )
    latex_lines += [r"\hline", r"\end{tabular}", r"\end{table}"]

    report = BenchmarkReport(
        results=obs_results,
        score=score,
        n_passed=n_passed,
        n_total=n_total,
        gamma_jet=system.gamma_jet,
        sigma_phi_satisfied=system._calibration.get("sigma_phi_frame_principle_satisfied", False),
        latex_table="\n".join(latex_lines),
        zenodo_json={
            "benchmark_score": score,
            "n_passed": n_passed,
            "n_total": n_total,
            "gamma_jet": system.gamma_jet,
            "observables": [
                {
                    "name": r.name,
                    "measured": r.measured,
                    "target": r.target,
                    "deviation_pct": r.deviation_pct,
                    "passed": r.passed,
                }
                for r in obs_results
            ],
        },
    )

    # Return as a plain dict for Diamond-Template compatibility
    return {
        "score": report.score,
        "n_passed": report.n_passed,
        "n_total": report.n_total,
        "passed": report.n_passed == report.n_total,
        "gamma_jet": report.gamma_jet,
        "sigma_phi_satisfied": report.sigma_phi_satisfied,
        "latex_table": report.latex_table,
        "zenodo_json": report.zenodo_json,
        "report_str": str(report),
        "observables": {
            r.name: {
                "measured": r.measured,
                "target": r.target,
                "tolerance": r.tolerance,
                "deviation_pct": r.deviation_pct,
                "passed": r.passed,
            }
            for r in obs_results
        },
    }
