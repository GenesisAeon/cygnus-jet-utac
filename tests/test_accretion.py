"""Tests for the accretion disk model."""

import pytest

from cygnus_jet_utac.accretion import AccretionDisk


@pytest.fixture
def disk() -> AccretionDisk:
    return AccretionDisk(seed=42)


class TestEddingtonLimits:
    def test_eddington_luminosity_positive(self, disk: AccretionDisk) -> None:
        assert disk.eddington_luminosity() > 0.0

    def test_eddington_luminosity_scales_with_mass(self) -> None:
        d1 = AccretionDisk(bh_mass=10.0 * 1.989e30)
        d2 = AccretionDisk(bh_mass=20.0 * 1.989e30)
        assert abs(d2.eddington_luminosity() / d1.eddington_luminosity() - 2.0) < 1e-9

    def test_eddington_accretion_rate_positive(self, disk: AccretionDisk) -> None:
        assert disk.eddington_accretion_rate() > 0.0

    def test_eddington_luminosity_for_21_msun(self, disk: AccretionDisk) -> None:
        """L_Edd = 4πGMm_p c/σ_T ≈ 2.64e32 W for 21 M☉ (~686 000 L☉)."""
        L = disk.eddington_luminosity()
        assert 1e32 < L < 1e33


class TestDiskLuminosity:
    def test_disk_luminosity_positive(self, disk: AccretionDisk) -> None:
        mdot = disk.eddington_accretion_rate() * 0.1
        assert disk.disk_luminosity(mdot) > 0.0

    def test_disk_luminosity_linear_in_mdot(self, disk: AccretionDisk) -> None:
        mdot = disk.eddington_accretion_rate() * 0.1
        L1 = disk.disk_luminosity(mdot)
        L2 = disk.disk_luminosity(2.0 * mdot)
        assert abs(L2 / L1 - 2.0) < 1e-9


class TestJetPower:
    def test_jet_power_10_percent_of_disk(self, disk: AccretionDisk) -> None:
        mdot = disk.eddington_accretion_rate() * 0.1
        P_jet = disk.jet_power(mdot, 0.10)
        P_disk = disk.disk_luminosity(mdot)
        assert abs(P_jet / P_disk - 1.0) < 1e-9  # both use η=0.10

    def test_jet_power_positive(self, disk: AccretionDisk) -> None:
        mdot = disk.eddington_accretion_rate() * 0.08
        assert disk.jet_power(mdot) > 0.0


class TestUTACNormalisation:
    def test_to_utac_H_at_eddington(self, disk: AccretionDisk) -> None:
        mdot_edd = disk.eddington_accretion_rate()
        H = disk.to_utac_H(mdot_edd)
        assert pytest.approx(1.0, abs=1e-9) == H

    def test_to_utac_H_zero_at_zero_mdot(self, disk: AccretionDisk) -> None:
        H = disk.to_utac_H(0.0)
        assert pytest.approx(0.0, abs=1e-12) == H

    def test_to_utac_H_clamped(self, disk: AccretionDisk) -> None:
        H = disk.to_utac_H(disk.eddington_accretion_rate() * 100.0)
        assert pytest.approx(1.0, abs=1e-9) == H

    def test_roundtrip(self, disk: AccretionDisk) -> None:
        mdot = disk.eddington_accretion_rate() * 0.08
        H = disk.to_utac_H(mdot)
        mdot_back = disk.from_utac_H(H)
        assert abs(mdot_back - mdot) / mdot < 1e-9

    def test_accretion_rate_modulated(self, disk: AccretionDisk) -> None:
        """accretion_rate should vary with time (orbital modulation)."""
        rates = [disk.accretion_rate(t) for t in range(0, int(disk.orbital_period), 3600)]
        assert max(rates) > min(rates)


class TestGeometry:
    def test_schwarzschild_radius_positive(self, disk: AccretionDisk) -> None:
        assert disk.schwarzschild_radius() > 0.0

    def test_isco_three_times_rs(self, disk: AccretionDisk) -> None:
        rs = disk.schwarzschild_radius()
        isco = disk.innermost_stable_orbit()
        assert abs(isco / rs - 3.0) < 1e-9

    def test_to_dict_keys(self, disk: AccretionDisk) -> None:
        d = disk.to_dict()
        assert "eddington_luminosity_W" in d
        assert "schwarzschild_radius_m" in d
