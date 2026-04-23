"""Tests for the relativistic jet propagation model."""

import math
import pytest
import numpy as np

from cygnus_jet_utac.jet import RelJet
from cygnus_jet_utac.orbital import CygnusOrbit
from cygnus_jet_utac.stellar_wind import StellarWindModel
from cygnus_jet_utac.constants import C_LIGHT, CYG_JET_VELOCITY, LIGHT_YEAR


@pytest.fixture
def jet() -> RelJet:
    return RelJet(seed=42)


@pytest.fixture
def orbit() -> CygnusOrbit:
    return CygnusOrbit()


@pytest.fixture
def wind() -> StellarWindModel:
    return StellarWindModel()


class TestKinematics:
    def test_beta_is_0p5(self, jet: RelJet) -> None:
        assert jet.beta == pytest.approx(0.5, abs=1e-9)

    def test_lorentz_factor(self, jet: RelJet) -> None:
        """γ = 1/√(1−0.25) = 1/√0.75 ≈ 1.1547."""
        expected = 1.0 / math.sqrt(1.0 - 0.5**2)
        assert jet.lorentz_factor() == pytest.approx(expected, rel=1e-6)

    def test_velocity_is_half_c(self, jet: RelJet) -> None:
        assert abs(jet.velocity_m_s - 0.5 * C_LIGHT) / C_LIGHT < 1e-9

    def test_lorentz_greater_than_one(self, jet: RelJet) -> None:
        assert jet.lorentz_factor() > 1.0


class TestMomentumFlux:
    def test_momentum_flux_positive(self, jet: RelJet) -> None:
        assert jet.jet_momentum_flux() > 0.0

    def test_momentum_flux_formula(self, jet: RelJet) -> None:
        expected = jet.jet_power / jet.velocity_m_s
        assert abs(jet.jet_momentum_flux() - expected) / expected < 1e-9


class TestPropagation:
    def test_initial_direction_along_z(self, jet: RelJet) -> None:
        d = jet.direction
        assert d[2] == pytest.approx(1.0, abs=1e-9)
        assert abs(d[0]) < 1e-12
        assert abs(d[1]) < 1e-12

    def test_propagation_increases_z(self, jet: RelJet, wind: StellarWindModel,
                                      orbit: CygnusOrbit) -> None:
        jet.reset()
        for i in range(10):
            jet.propagate(3600.0, wind, orbit, i * 3600.0)
        # Jet should have moved in z-direction predominantly
        pos = jet.position
        assert pos[2] > 0.0

    def test_direction_is_unit_vector(self, jet: RelJet, wind: StellarWindModel,
                                       orbit: CygnusOrbit) -> None:
        jet.reset()
        for i in range(100):
            jet.propagate(3600.0, wind, orbit, i * 3600.0)
        d = jet.direction
        assert abs(np.linalg.norm(d) - 1.0) < 1e-9

    def test_extent_increases_monotonically(self, jet: RelJet, wind: StellarWindModel,
                                             orbit: CygnusOrbit) -> None:
        jet.reset()
        extents = []
        for i in range(20):
            jet.propagate(3600.0, wind, orbit, i * 3600.0)
            extents.append(jet.extent())
        assert extents == sorted(extents)

    def test_reset_clears_position(self, jet: RelJet) -> None:
        jet.reset()
        assert np.allclose(jet.position, np.zeros(3))
        assert jet.extent() == 0.0


class TestObservables:
    def test_extent_ly_positive_after_propagation(self, jet: RelJet,
                                                    wind: StellarWindModel,
                                                    orbit: CygnusOrbit) -> None:
        jet.reset()
        for i in range(100):
            jet.propagate(86400.0, wind, orbit, i * 86400.0)
        assert jet.extent_ly() > 0.0

    def test_to_dict_keys(self, jet: RelJet) -> None:
        d = jet.to_dict()
        assert "beta" in d
        assert "gamma_lorentz" in d
        assert "jet_power_W" in d
        assert "extent_ly" in d
