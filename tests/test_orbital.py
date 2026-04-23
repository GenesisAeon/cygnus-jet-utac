"""Tests for the Cygnus X-1 binary orbital dynamics."""

import math
import pytest
import numpy as np

from cygnus_jet_utac.orbital import CygnusOrbit
from cygnus_jet_utac.constants import CYG_ORBITAL_PERIOD, CYG_BH_MASS, CYG_COMPANION_MASS


@pytest.fixture
def orbit() -> CygnusOrbit:
    return CygnusOrbit()


class TestOrbitalPeriod:
    def test_period_is_5p6_days(self, orbit: CygnusOrbit) -> None:
        assert abs(orbit.period / 86400.0 - 5.6) < 0.001

    def test_phase_zero_at_t0(self, orbit: CygnusOrbit) -> None:
        assert orbit.phase(0.0) == pytest.approx(0.0, abs=1e-12)

    def test_phase_2pi_after_one_period(self, orbit: CygnusOrbit) -> None:
        # After one period, phase should return to 0 (mod 2π)
        phi = orbit.phase(CYG_ORBITAL_PERIOD)
        assert phi == pytest.approx(0.0, abs=1e-9)

    def test_phase_pi_at_half_period(self, orbit: CygnusOrbit) -> None:
        phi = orbit.phase(0.5 * CYG_ORBITAL_PERIOD)
        assert phi == pytest.approx(math.pi, abs=1e-9)

    def test_phase_range(self, orbit: CygnusOrbit) -> None:
        for t in np.linspace(0, 10 * CYG_ORBITAL_PERIOD, 1000):
            phi = orbit.phase(t)
            assert 0.0 <= phi <= 2.0 * math.pi + 1e-12


class TestSeparation:
    def test_separation_positive(self, orbit: CygnusOrbit) -> None:
        assert orbit.separation(0.0) > 0.0

    def test_separation_constant(self, orbit: CygnusOrbit) -> None:
        """Circular orbit has constant separation."""
        seps = [orbit.separation(t) for t in [0, 1e5, 1e6, 1e7]]
        assert all(abs(s - seps[0]) < 1.0 for s in seps)

    def test_separation_in_au_range(self, orbit: CygnusOrbit) -> None:
        """Cygnus X-1 separation ~0.2 AU."""
        sep_AU = orbit.separation(0.0) / 1.496e11
        assert 0.05 < sep_AU < 1.0


class TestPositions:
    def test_bh_position_at_t0(self, orbit: CygnusOrbit) -> None:
        pos = orbit.position_bh(0.0)
        assert pos.shape == (3,)
        assert pos[2] == pytest.approx(0.0, abs=1e-10)

    def test_companion_opposite_bh(self, orbit: CygnusOrbit) -> None:
        """BH and companion should be on opposite sides of CoM."""
        pos_bh = orbit.position_bh(0.0)
        pos_comp = orbit.position_companion(0.0)
        # Their x-components should have opposite signs
        assert pos_bh[0] * pos_comp[0] <= 0.0

    def test_positions_sum_to_zero_weighted(self, orbit: CygnusOrbit) -> None:
        """M1·r1 + M2·r2 = 0 (centre of mass) to floating-point precision."""
        for t in [0.0, 1e5, 3e5]:
            p1 = orbit.position_bh(t)
            p2 = orbit.position_companion(t)
            com = CYG_BH_MASS * p1 + CYG_COMPANION_MASS * p2
            # Scale by M*r ~ 1e42; floating-point residual ~ 1e42 * eps ~ 1e26
            scale = np.linalg.norm(CYG_BH_MASS * p1)
            relative_error = np.linalg.norm(com) / scale if scale > 0 else 0.0
            assert relative_error < 1e-10


class TestWindIncidentAngle:
    def test_angle_in_range(self, orbit: CygnusOrbit) -> None:
        for t in np.linspace(0, CYG_ORBITAL_PERIOD, 100):
            theta = orbit.wind_incident_angle(t)
            assert 0.0 <= theta <= math.pi + 1e-9

    def test_angle_periodic(self, orbit: CygnusOrbit) -> None:
        t = 1.234e6
        theta1 = orbit.wind_incident_angle(t)
        theta2 = orbit.wind_incident_angle(t + CYG_ORBITAL_PERIOD)
        assert abs(theta1 - theta2) < 1e-9


class TestMiscellaneous:
    def test_angular_velocity_positive(self, orbit: CygnusOrbit) -> None:
        assert orbit.angular_velocity > 0.0

    def test_orbital_velocity_positive(self, orbit: CygnusOrbit) -> None:
        assert orbit.orbital_velocity() > 0.0

    def test_roche_lobe_positive(self, orbit: CygnusOrbit) -> None:
        assert orbit.roche_lobe_radius() > 0.0

    def test_to_dict_keys(self, orbit: CygnusOrbit) -> None:
        d = orbit.to_dict()
        assert "period_days" in d
        assert "separation_AU" in d
