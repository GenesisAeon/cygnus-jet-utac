"""Tests for the benchmark module validating against Prabu et al. 2026."""

import pytest
from cygnus_jet_utac.benchmark import (
    PRABU_2026_TARGETS,
    run_benchmark,
    BenchmarkReport,
    ObservableResult,
)
from cygnus_jet_utac.system import CygnusJetUTAC


@pytest.fixture(scope="module")
def benchmarked_system() -> CygnusJetUTAC:
    s = CygnusJetUTAC(dt=21600.0, seed=42)
    s.run_cycle(duration_years=18.0)
    return s


class TestPrabu2026Targets:
    def test_targets_dict_has_required_keys(self) -> None:
        required = {
            "jet_power_W", "jet_velocity_c", "accretion_efficiency",
            "jet_extent_ly", "orbital_period_days", "dance_events_per_year",
        }
        assert required.issubset(PRABU_2026_TARGETS.keys())

    def test_targets_have_correct_values(self) -> None:
        assert PRABU_2026_TARGETS["jet_power_W"][0] == pytest.approx(3.846e30, rel=1e-3)
        assert PRABU_2026_TARGETS["jet_velocity_c"][0] == pytest.approx(0.50, rel=1e-6)
        assert PRABU_2026_TARGETS["accretion_efficiency"][0] == pytest.approx(0.10, rel=1e-6)
        assert PRABU_2026_TARGETS["jet_extent_ly"][0] == pytest.approx(16.0, rel=1e-6)
        assert PRABU_2026_TARGETS["orbital_period_days"][0] == pytest.approx(5.6, rel=1e-6)

    def test_tolerances_positive(self) -> None:
        for name, (val, tol) in PRABU_2026_TARGETS.items():
            assert tol > 0.0, f"Tolerance for {name} must be positive"


class TestRunBenchmark:
    def test_returns_dict(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert isinstance(result, dict)

    def test_has_score(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "score" in result
        assert 0.0 <= result["score"] <= 100.0

    def test_has_n_passed(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "n_passed" in result
        assert 0 <= result["n_passed"] <= result["n_total"]

    def test_has_observables(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "observables" in result
        assert len(result["observables"]) == len(PRABU_2026_TARGETS)

    def test_has_gamma_jet(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "gamma_jet" in result
        assert abs(result["gamma_jet"] - 0.0456) < 0.001

    def test_has_latex_table(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "latex_table" in result
        assert r"\begin{table}" in result["latex_table"]

    def test_has_zenodo_json(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "zenodo_json" in result
        assert "benchmark_score" in result["zenodo_json"]

    def test_has_report_str(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        assert "report_str" in result
        assert "Prabu" in result["report_str"]

    def test_jet_velocity_within_tolerance(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        obs = result["observables"]["jet_velocity_c"]
        assert obs["passed"], (
            f"jet_velocity_c deviation {obs['deviation_pct']:.1f}% > "
            f"{obs['tolerance']*100:.0f}% tolerance"
        )

    def test_orbital_period_within_tolerance(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        obs = result["observables"]["orbital_period_days"]
        assert obs["passed"], (
            f"orbital_period deviation {obs['deviation_pct']:.1f}% > "
            f"{obs['tolerance']*100:.0f}% tolerance"
        )

    def test_accretion_efficiency_within_tolerance(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        obs = result["observables"]["accretion_efficiency"]
        assert obs["passed"], (
            f"efficiency deviation {obs['deviation_pct']:.1f}% > "
            f"{obs['tolerance']*100:.0f}% tolerance"
        )

    def test_jet_power_within_tolerance(self, benchmarked_system: CygnusJetUTAC) -> None:
        result = run_benchmark(benchmarked_system)
        obs = result["observables"]["jet_power_W"]
        assert obs["passed"], (
            f"jet_power_W deviation {obs['deviation_pct']:.1f}% > "
            f"{obs['tolerance']*100:.0f}% tolerance"
        )

    def test_score_above_threshold(self, benchmarked_system: CygnusJetUTAC) -> None:
        """Overall benchmark score must reach ≥ 80/100."""
        result = run_benchmark(benchmarked_system)
        assert result["score"] >= 80.0, (
            f"Benchmark score {result['score']:.1f} < 80.  Report:\n"
            f"{result['report_str']}"
        )

    def test_auto_runs_if_no_results(self) -> None:
        """run_benchmark should auto-call run_cycle if no results present."""
        s = CygnusJetUTAC(dt=86400.0, seed=99)
        result = run_benchmark(s)
        assert "score" in result
