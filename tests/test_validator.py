"""Tests for calibration and parameter validation."""

import math

from cygnus_jet_utac.constants import (
    CYG_ACCRETION_EFFICIENCY,
    UTAC_R_DEFAULT,
    UTAC_SIGMA_DEFAULT,
)
from cygnus_jet_utac.efficiency import calibrate_gamma_jet


def test_calibrate_returns_dict() -> None:
    result = calibrate_gamma_jet()
    assert isinstance(result, dict)


def test_calibrate_gamma_jet_value() -> None:
    result = calibrate_gamma_jet(eta=0.10, sigma=2.2)
    expected = math.atanh(0.10) / 2.2
    assert abs(result["gamma_jet"] - expected) < 1e-6


def test_frame_principle_satisfied() -> None:
    result = calibrate_gamma_jet(eta=0.10, sigma=2.2)
    assert result["sigma_phi_frame_principle_satisfied"] is True


def test_sigma_phi_ratio_in_range() -> None:
    result = calibrate_gamma_jet(eta=0.10, sigma=2.2)
    ratio = result["sigma_phi_ratio"]
    assert 0.0 < ratio <= 1.0


def test_utac_r_default_positive() -> None:
    assert UTAC_R_DEFAULT > 0.0


def test_utac_sigma_default_positive() -> None:
    assert UTAC_SIGMA_DEFAULT > 0.0


def test_accretion_efficiency_in_range() -> None:
    assert 0.0 < CYG_ACCRETION_EFFICIENCY < 1.0
