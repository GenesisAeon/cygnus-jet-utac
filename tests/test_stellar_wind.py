"""Tests for the Blue Supergiant stellar wind model."""

import math
import pytest
import numpy as np

from cygnus_jet_utac.stellar_wind import StellarWindModel
from cygnus_jet_utac.constants import CYG_ORBITAL_SEPARATION, CYG_COMPANION_RADIUS


@pytest.fixture
def wind() -> StellarWindModel:
    return StellarWindModel()


class TestVelocityLaw:
    def test_velocity_zero_at_surface(self, wind: StellarWindModel) -> None:
        v = wind.velocity(CYG_COMPANION_RADIUS)
        assert v == pytest.approx(0.0, abs=1e-3)

    def test_velocity_increases_with_r(self, wind: StellarWindModel) -> None:
        # Use multiples of R_star, sorted in ascending order
        r_values = sorted([
            1.1 * CYG_COMPANION_RADIUS, 2.0 * CYG_COMPANION_RADIUS,
            5.0 * CYG_COMPANION_RADIUS, 10.0 * CYG_COMPANION_RADIUS,
        ])
        vels = [wind.velocity(r) for r in r_values]
        assert vels == sorted(vels)

    def test_velocity_approaches_vinf(self, wind: StellarWindModel) -> None:
        """At very large r, v → v_inf."""
        v_large = wind.velocity(1e6 * CYG_COMPANION_RADIUS)
        assert abs(v_large / wind.v_inf - 1.0) < 0.001

    def test_velocity_positive_outside_star(self, wind: StellarWindModel) -> None:
        r = 2.0 * CYG_COMPANION_RADIUS
        assert wind.velocity(r) > 0.0

    def test_velocity_zero_inside_star(self, wind: StellarWindModel) -> None:
        v = wind.velocity(0.5 * CYG_COMPANION_RADIUS)
        assert v == 0.0


class TestDensity:
    def test_density_decreases_with_r(self, wind: StellarWindModel) -> None:
        r_values = sorted([
            1.5 * CYG_COMPANION_RADIUS, 3.0 * CYG_COMPANION_RADIUS,
            8.0 * CYG_COMPANION_RADIUS, 15.0 * CYG_COMPANION_RADIUS,
        ])
        densities = [wind.mass_density(r) for r in r_values]
        assert densities == sorted(densities, reverse=True)

    def test_density_positive(self, wind: StellarWindModel) -> None:
        assert wind.mass_density(CYG_ORBITAL_SEPARATION) > 0.0

    def test_number_density_positive(self, wind: StellarWindModel) -> None:
        assert wind.number_density(CYG_ORBITAL_SEPARATION) > 0.0


class TestRamPressure:
    def test_ram_pressure_positive(self, wind: StellarWindModel) -> None:
        p = wind.ram_pressure(CYG_ORBITAL_SEPARATION)
        assert p > 0.0

    def test_ram_pressure_decreases_with_r(self, wind: StellarWindModel) -> None:
        p1 = wind.ram_pressure(1.5 * CYG_COMPANION_RADIUS)
        p2 = wind.ram_pressure(CYG_ORBITAL_SEPARATION)
        assert p1 > p2


class TestCREPComponent:
    def test_crep_R_in_range(self, wind: StellarWindModel) -> None:
        jet_pos = np.array([0.0, 0.0, 1e12])
        for t in np.linspace(0, wind.orbital_period * 3, 50):
            R = wind.crep_R_component(t, jet_pos)
            assert 0.0 <= R <= 1.0, f"R={R} out of [0,1] at t={t}"

    def test_crep_R_periodic(self, wind: StellarWindModel) -> None:
        """R should be approximately periodic with orbital period."""
        jet_pos = np.array([0.0, 0.0, 1e12])
        T = wind.orbital_period
        R0 = wind.crep_R_component(0.0, jet_pos)
        R1 = wind.crep_R_component(T, jet_pos)
        assert abs(R0 - R1) < 0.01  # within 1% after one period


class TestWindForce:
    def test_wind_force_shape(self, wind: StellarWindModel) -> None:
        d = np.array([0.0, 0.0, 1.0])
        f = wind.wind_force_on_jet(CYG_ORBITAL_SEPARATION, d, math.pi / 4)
        assert f.shape == (3,)

    def test_wind_force_zero_along_jet(self, wind: StellarWindModel) -> None:
        """Force perpendicular to jet; component along jet must be zero."""
        d = np.array([0.0, 0.0, 1.0])
        f = wind.wind_force_on_jet(CYG_ORBITAL_SEPARATION, d, math.pi / 4)
        assert abs(np.dot(f, d)) < 1e-10


class TestToDict:
    def test_to_dict_keys(self, wind: StellarWindModel) -> None:
        d = wind.to_dict()
        assert "v_inf_m_s" in d
        assert "beta" in d
        assert "mdot_Msun_yr" in d
