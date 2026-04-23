"""Mandatory Diamond-Template contract tests for CygnusJetUTAC.

These tests are non-negotiable for GenesisAeon Package 17 certification.
Every method in the Diamond-Template interface must exist, return the correct
type, and carry the required keys.
"""

import pytest

from cygnus_jet_utac.system import CygnusJetUTAC


@pytest.fixture(scope="module")
def system() -> CygnusJetUTAC:
    """Short simulation (1 year, 6-hour steps) for fast tests."""
    s = CygnusJetUTAC(dt=21600.0, seed=42)
    s.run_cycle(duration_years=1.0)
    return s


@pytest.fixture(scope="module")
def fresh_system() -> CygnusJetUTAC:
    """Fresh system (no run_cycle) for state-inspection tests."""
    return CygnusJetUTAC(dt=86400.0, seed=0)


# ── run_cycle ─────────────────────────────────────────────────────────────────

def test_run_cycle_returns_dict(system: CygnusJetUTAC) -> None:
    result = system.run_cycle(duration_years=0.5)
    assert isinstance(result, dict)


def test_run_cycle_has_jet_power(system: CygnusJetUTAC) -> None:
    result = system.run_cycle(duration_years=0.5)
    assert "jet_power_W" in result
    assert result["jet_power_W"] > 0.0


def test_run_cycle_has_gamma_jet(system: CygnusJetUTAC) -> None:
    result = system.run_cycle(duration_years=0.5)
    assert "gamma_jet" in result
    assert 0.0 < result["gamma_jet"] < 1.0


def test_run_cycle_has_phase_events_list(system: CygnusJetUTAC) -> None:
    result = system.run_cycle(duration_years=0.5)
    assert isinstance(result["phase_events"], list)


def test_run_cycle_reproducible(fresh_system: CygnusJetUTAC) -> None:
    s1 = CygnusJetUTAC(dt=86400.0, seed=42)
    s2 = CygnusJetUTAC(dt=86400.0, seed=42)
    r1 = s1.run_cycle(duration_years=0.25)
    r2 = s2.run_cycle(duration_years=0.25)
    assert abs(r1["jet_power_W"] - r2["jet_power_W"]) < 1.0


# ── get_crep_state ────────────────────────────────────────────────────────────

def test_get_crep_state_returns_dict(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_crep_state()
    assert isinstance(state, dict)


def test_get_crep_state_keys(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_crep_state()
    for key in ("C", "R", "E", "P", "Gamma"):
        assert key in state, f"Missing CREP key: {key}"


def test_get_crep_state_values_in_range(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_crep_state()
    for key in ("C", "R", "E", "P", "Gamma"):
        assert 0.0 <= state[key] <= 1.0, f"CREP {key} = {state[key]} out of [0,1]"


# ── get_utac_state ────────────────────────────────────────────────────────────

def test_get_utac_state_returns_dict(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_utac_state()
    assert isinstance(state, dict)


def test_get_utac_state_keys(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_utac_state()
    for key in ("H", "dH_dt", "H_star", "K_eff"):
        assert key in state, f"Missing UTAC key: {key}"


def test_get_utac_state_H_nonnegative(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_utac_state()
    assert state["H"] >= 0.0


def test_get_utac_state_H_star_positive(fresh_system: CygnusJetUTAC) -> None:
    state = fresh_system.get_utac_state()
    assert state["H_star"] >= 0.0


# ── get_phase_events ──────────────────────────────────────────────────────────

def test_get_phase_events_is_list(system: CygnusJetUTAC) -> None:
    events = system.get_phase_events()
    assert isinstance(events, list)


def test_get_phase_events_structure(system: CygnusJetUTAC) -> None:
    # Run longer to ensure at least one event
    s = CygnusJetUTAC(dt=21600.0, seed=42)
    s.run_cycle(duration_years=3.0)
    events = s.get_phase_events()
    if events:
        ev = events[0]
        assert "t" in ev
        assert "deflection_angle_deg" in ev
        assert "crep_gamma" in ev


# ── to_zenodo_record ──────────────────────────────────────────────────────────

def test_to_zenodo_record_returns_dict(fresh_system: CygnusJetUTAC) -> None:
    record = fresh_system.to_zenodo_record()
    assert isinstance(record, dict)


def test_to_zenodo_record_has_doi(fresh_system: CygnusJetUTAC) -> None:
    record = fresh_system.to_zenodo_record()
    assert "doi" in record
    assert record["doi"] == "10.5281/zenodo.19645351"


def test_to_zenodo_record_has_reference_doi(fresh_system: CygnusJetUTAC) -> None:
    record = fresh_system.to_zenodo_record()
    assert "reference_doi" in record
    assert "10.1038" in record["reference_doi"]


def test_to_zenodo_record_has_parameters(fresh_system: CygnusJetUTAC) -> None:
    record = fresh_system.to_zenodo_record()
    assert "parameters" in record
    assert "bh_mass_Msun" in record["parameters"]


def test_to_zenodo_record_has_calibration(fresh_system: CygnusJetUTAC) -> None:
    record = fresh_system.to_zenodo_record()
    assert "calibration" in record
    assert "gamma_jet" in record["calibration"]


# ── __repr__ ──────────────────────────────────────────────────────────────────

def test_repr_contains_gamma_jet(fresh_system: CygnusJetUTAC) -> None:
    r = repr(fresh_system)
    assert "Γ_jet" in r
    assert "CygnusJetUTAC" in r
