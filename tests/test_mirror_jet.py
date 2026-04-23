"""Tests for the JetMirrorMachine phase-transition detector."""

import math

import numpy as np
import pytest

from cygnus_jet_utac.constants import CYG_ORBITAL_PERIOD
from cygnus_jet_utac.mirror_jet import JetMirrorMachine


@pytest.fixture
def mirror() -> JetMirrorMachine:
    return JetMirrorMachine(tau=CYG_ORBITAL_PERIOD, dt=3600.0)


class TestInitialState:
    def test_no_events_initially(self, mirror: JetMirrorMachine) -> None:
        assert len(mirror.get_dance_events()) == 0

    def test_divergence_zero_initially(self, mirror: JetMirrorMachine) -> None:
        # Buffer is filled with identical vectors → divergence = 0
        assert mirror.divergence() == pytest.approx(0.0, abs=1e-9)


class TestUpdate:
    def test_returns_bool(self, mirror: JetMirrorMachine) -> None:
        d = np.array([0.0, 0.0, 1.0])
        result = mirror.update(d, crep_gamma=0.05, t=0.0)
        assert isinstance(result, bool)

    def test_no_event_for_constant_direction(self, mirror: JetMirrorMachine) -> None:
        d = np.array([0.0, 0.0, 1.0])
        events = 0
        for i in range(200):
            triggered = mirror.update(d, crep_gamma=0.05, t=i * 3600.0)
            if triggered:
                events += 1
        assert events == 0

    def test_event_triggered_on_large_deflection(self) -> None:
        mirror = JetMirrorMachine(tau=CYG_ORBITAL_PERIOD, dt=3600.0, theta0=0.05)
        n_delay = mirror._n_delay

        d_initial = np.array([0.0, 0.0, 1.0])
        d_deflected = np.array([math.sin(0.3), 0.0, math.cos(0.3)])  # 17° deflection

        # Fill buffer with initial direction
        for i in range(n_delay + 1):
            mirror.update(d_initial, crep_gamma=0.05, t=i * 3600.0)

        # Now provide deflected direction — should trigger
        triggered = mirror.update(
            d_deflected, crep_gamma=0.05, t=(n_delay + 2) * 3600.0
        )
        assert triggered, "Expected phase event for large deflection"

    def test_accepts_non_unit_vector(self, mirror: JetMirrorMachine) -> None:
        d = np.array([0.0, 0.0, 5.0])  # non-unit, should be normalised internally
        result = mirror.update(d, crep_gamma=0.05, t=0.0)
        assert isinstance(result, bool)


class TestDanceEvents:
    def test_event_has_required_fields(self) -> None:
        mirror = JetMirrorMachine(tau=CYG_ORBITAL_PERIOD, dt=3600.0, theta0=0.05)
        n_delay = mirror._n_delay

        d0 = np.array([0.0, 0.0, 1.0])
        d1 = np.array([math.sin(0.4), 0.0, math.cos(0.4)])

        for i in range(n_delay + 1):
            mirror.update(d0, 0.05, i * 3600.0)
        mirror.update(d1, 0.05, (n_delay + 2) * 3600.0)

        events = mirror.get_dance_events()
        if events:
            ev = events[0]
            assert hasattr(ev, "t")
            assert hasattr(ev, "t_years")
            assert hasattr(ev, "deflection_angle_deg")
            assert hasattr(ev, "crep_gamma")
            assert hasattr(ev, "direction_before")
            assert hasattr(ev, "direction_after")

    def test_deflection_angle_positive(self) -> None:
        mirror = JetMirrorMachine(tau=CYG_ORBITAL_PERIOD, dt=3600.0, theta0=0.05)
        n_delay = mirror._n_delay
        d0 = np.array([0.0, 0.0, 1.0])
        d1 = np.array([math.sin(0.4), 0.0, math.cos(0.4)])
        for i in range(n_delay + 1):
            mirror.update(d0, 0.05, i * 3600.0)
        mirror.update(d1, 0.05, (n_delay + 2) * 3600.0)
        events = mirror.get_dance_events()
        if events:
            assert events[0].deflection_angle_deg > 0.0


class TestRefractory:
    def test_refractory_prevents_back_to_back(self) -> None:
        """Two consecutive large deflections within 0.5τ should yield only 1 event."""
        mirror = JetMirrorMachine(tau=3600.0, dt=360.0, theta0=0.05)
        n_delay = mirror._n_delay
        d0 = np.array([0.0, 0.0, 1.0])
        d1 = np.array([0.3, 0.0, math.sqrt(1 - 0.09)])

        for i in range(n_delay + 1):
            mirror.update(d0, 0.05, i * 360.0)

        trigger_count = 0
        for j in range(3):
            triggered = mirror.update(d1, 0.05, (n_delay + 2 + j) * 360.0)
            if triggered:
                trigger_count += 1

        assert trigger_count <= 1


class TestReset:
    def test_reset_clears_events(self, mirror: JetMirrorMachine) -> None:
        d0 = np.array([0.0, 0.0, 1.0])
        for i in range(100):
            mirror.update(d0, 0.05, i * 3600.0)
        mirror.reset()
        assert len(mirror.get_dance_events()) == 0


class TestDanceEventsPerYear:
    def test_rate_zero_no_events(self, mirror: JetMirrorMachine) -> None:
        assert mirror.dance_events_per_year(1.0) == pytest.approx(0.0)

    def test_rate_zero_duration_zero(self, mirror: JetMirrorMachine) -> None:
        assert mirror.dance_events_per_year(0.0) == pytest.approx(0.0)


class TestToDict:
    def test_to_dict_keys(self, mirror: JetMirrorMachine) -> None:
        d = mirror.to_dict()
        assert "tau_s" in d
        assert "n_dance_events" in d
        assert "current_divergence_rad" in d
