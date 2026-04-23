"""Tests for the Γ_jet calibration module — the central scientific result."""

import math

import numpy as np
import pytest

from cygnus_jet_utac.efficiency import (
    calibrate_gamma_jet,
    efficiency_from_gamma,
    gamma_scan,
)


class TestCalibrateGammaJet:
    def test_gamma_jet_from_eta_10_percent(self) -> None:
        """Γ_jet = arctanh(0.10) / 2.2 ≈ 0.0456."""
        result = calibrate_gamma_jet(eta=0.10)
        assert abs(result["gamma_jet"] - 0.0456) < 0.001

    def test_gamma_jet_exact_formula(self) -> None:
        """Verify Γ_jet = arctanh(η) / σ exactly."""
        eta, sigma = 0.10, 2.2
        expected = math.atanh(eta) / sigma
        result = calibrate_gamma_jet(eta=eta, sigma=sigma)
        assert abs(result["gamma_jet"] - expected) < 1e-10

    def test_utac_fixed_point_check(self) -> None:
        """K·tanh(σ·Γ_jet) must equal η·K."""
        result = calibrate_gamma_jet(eta=0.10, K=1.0)
        assert abs(result["utac_fixed_point_check"] - 0.10) < 1e-10

    def test_efficiency_check_roundtrip(self) -> None:
        """Recovered η must equal input η to floating-point precision."""
        for eta in [0.05, 0.10, 0.20, 0.50]:
            result = calibrate_gamma_jet(eta=eta)
            assert abs(result["efficiency_check"] - eta) < 1e-10, (
                f"Roundtrip failed for η={eta}: got {result['efficiency_check']}"
            )

    def test_sigma_phi_near_1_over_16(self) -> None:
        """σ_Φ,min should be ≤ 1/16 for Frame Principle to hold."""
        result = calibrate_gamma_jet(eta=0.10)
        assert result["sigma_phi_ratio"] <= 1.0, (
            f"Frame Principle violated: ratio = {result['sigma_phi_ratio']:.4f} > 1"
        )

    def test_sigma_phi_ratio_reasonable(self) -> None:
        """σ_Φ ratio should be in a physically plausible range (0.5 – 1.0)."""
        result = calibrate_gamma_jet(eta=0.10)
        assert 0.3 < result["sigma_phi_ratio"] < 1.1

    def test_frame_principle_flag(self) -> None:
        """sigma_phi_frame_principle_satisfied must be True for η=10%."""
        result = calibrate_gamma_jet(eta=0.10)
        assert result["sigma_phi_frame_principle_satisfied"] is True

    def test_verbose_report_present(self) -> None:
        result = calibrate_gamma_jet(eta=0.10, verbose=True)
        assert "report" in result
        assert "Γ_jet" in result["report"]

    def test_no_report_without_verbose(self) -> None:
        result = calibrate_gamma_jet(eta=0.10, verbose=False)
        assert "report" not in result

    def test_invalid_eta_raises(self) -> None:
        with pytest.raises(ValueError):
            calibrate_gamma_jet(eta=0.0)
        with pytest.raises(ValueError):
            calibrate_gamma_jet(eta=1.0)
        with pytest.raises(ValueError):
            calibrate_gamma_jet(eta=-0.1)

    def test_interpretation_in_result(self) -> None:
        result = calibrate_gamma_jet(eta=0.10)
        assert "interpretation" in result
        assert len(result["interpretation"]) > 10

    def test_gamma_jet_subcritical(self) -> None:
        """Cygnus X-1 jet must be sub-critical (Γ < 0.5)."""
        result = calibrate_gamma_jet(eta=0.10)
        assert result["gamma_jet"] < 0.5


class TestEfficiencyFromGamma:
    def test_roundtrip_with_calibration(self) -> None:
        """efficiency_from_gamma(gamma_jet) ≈ 0.10."""
        from cygnus_jet_utac.efficiency import calibrate_gamma_jet
        gamma_jet = calibrate_gamma_jet(eta=0.10)["gamma_jet"]
        eta_back = efficiency_from_gamma(gamma_jet)
        assert abs(eta_back - 0.10) < 1e-9

    def test_monotone_increasing(self) -> None:
        """Higher Γ → higher efficiency (tanh is monotone)."""
        etas = [efficiency_from_gamma(g) for g in [0.01, 0.1, 0.5, 1.0, 2.0]]
        assert etas == sorted(etas)


class TestGammaScan:
    def test_returns_two_arrays(self) -> None:
        eta_arr, gamma_arr = gamma_scan()
        assert len(eta_arr) == len(gamma_arr)
        assert len(eta_arr) > 1

    def test_custom_range(self) -> None:
        eta_arr, gamma_arr = gamma_scan(np.array([0.05, 0.10, 0.20]))
        assert len(eta_arr) == 3

    def test_gamma_increases_with_eta(self) -> None:
        eta_arr, gamma_arr = gamma_scan()
        assert all(gamma_arr[i] < gamma_arr[i + 1] for i in range(len(gamma_arr) - 1))
