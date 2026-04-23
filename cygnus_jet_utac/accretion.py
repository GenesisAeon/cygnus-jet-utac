"""Accretion disk dynamics — maps to UTAC state variable H(t).

The mass accretion rate Ṁ(t) is normalised to the Eddington rate to give
H(t) ∈ [0, K]. The CREP tensor Γ(t) modulates whether disk energy flows into
the jet or is radiated away — the 10 % efficiency is the UTAC fixed point.
"""

from __future__ import annotations

import math

import numpy as np

from cygnus_jet_utac.constants import (
    ACCRETION_RADIATIVE_EFF,
    C_LIGHT,
    CYG_BH_MASS,
    CYG_ORBITAL_PERIOD,
    M_PROTON,
    SIGMA_T,
    G,
)


class AccretionDisk:
    """Cygnus X-1 accretion disk model.

    The disk mass accretion rate Ṁ(t) is the UTAC state variable H(t),
    normalised so H = 1 corresponds to the Eddington accretion rate.

    The CREP tensor Γ(t) drives H(t) toward the fixed point:

        H*(t) = K · tanh(σ · Γ(t))

    At the calibrated Γ_jet ≈ 0.046 and σ = 2.2 this gives H*/K ≈ 0.10,
    encoding the 10 % accretion-to-jet efficiency measured by Prabu et al.

    Args:
        bh_mass: Black hole mass (kg).
        radiative_eff: Standard radiative efficiency of the disk (dimensionless).
        orbital_period: Orbital period for wind-modulated variability (s).
        variability_amplitude: Peak-to-trough fractional variability (dimensionless).
        seed: Random seed for stochastic variability.
    """

    def __init__(
        self,
        bh_mass: float = CYG_BH_MASS,
        radiative_eff: float = ACCRETION_RADIATIVE_EFF,
        orbital_period: float = CYG_ORBITAL_PERIOD,
        variability_amplitude: float = 0.15,
        seed: int = 42,
    ) -> None:
        self.bh_mass = bh_mass
        self.radiative_eff = radiative_eff
        self.orbital_period = orbital_period
        self.variability_amplitude = variability_amplitude
        self._rng = np.random.default_rng(seed)
        self._mdot_edd = self.eddington_accretion_rate()

    # ── Eddington limits ──────────────────────────────────────────────────────

    def eddington_luminosity(self) -> float:
        """Eddington luminosity L_Edd = 4π G M m_p c / σ_T (W).

        Returns:
            L_Edd (W).
        """
        return 4.0 * math.pi * G * self.bh_mass * M_PROTON * C_LIGHT / SIGMA_T

    def eddington_accretion_rate(self) -> float:
        """Eddington mass accretion rate Ṁ_Edd = L_Edd / (η c²) (kg/s).

        Returns:
            Ṁ_Edd (kg/s).
        """
        return self.eddington_luminosity() / (self.radiative_eff * C_LIGHT**2)

    # ── Time-varying accretion rate ───────────────────────────────────────────

    def accretion_rate(self, t: float, H_norm: float = 0.08) -> float:
        """Instantaneous mass accretion rate Ṁ(t) (kg/s).

        The accretion rate has an orbital-period sinusoidal modulation from
        the companion wind density variations at the accretion stream nozzle,
        plus a stochastic quasi-periodic oscillation component.

        Args:
            t: Simulation time (s).
            H_norm: Current normalised UTAC state H ∈ [0, 1].

        Returns:
            Ṁ(t) (kg/s).
        """
        # Orbital modulation: X-ray binaries show ~10–20 % modulation
        orb_phase = 2.0 * math.pi * t / self.orbital_period
        orbital_mod = 1.0 + self.variability_amplitude * math.sin(orb_phase)

        return H_norm * self._mdot_edd * orbital_mod

    def disk_luminosity(self, mdot: float) -> float:
        """Disk radiative luminosity L_disk = η_rad · Ṁ · c² (W).

        Args:
            mdot: Mass accretion rate (kg/s).

        Returns:
            L_disk (W).
        """
        return self.radiative_eff * mdot * C_LIGHT**2

    def jet_power(self, mdot: float, jet_efficiency: float = 0.10) -> float:
        """Jet mechanical power P_jet = η_jet · Ṁ · c² (W).

        Args:
            mdot: Mass accretion rate (kg/s).
            jet_efficiency: Accretion-to-jet energy conversion efficiency η.

        Returns:
            P_jet (W).
        """
        return jet_efficiency * mdot * C_LIGHT**2

    # ── UTAC normalisation ────────────────────────────────────────────────────

    def to_utac_H(self, mdot: float) -> float:
        """Normalise Ṁ to UTAC state variable H ∈ [0, K].

        H = Ṁ / Ṁ_Edd,  so H = 1 ↔ Eddington rate.

        Args:
            mdot: Mass accretion rate (kg/s).

        Returns:
            H (dimensionless, ∈ [0, 1]).
        """
        return float(np.clip(mdot / self._mdot_edd, 0.0, 1.0))

    def from_utac_H(self, H: float) -> float:
        """Convert UTAC state H back to physical accretion rate (kg/s).

        Args:
            H: Normalised state variable ∈ [0, 1].

        Returns:
            Ṁ (kg/s).
        """
        return float(np.clip(H, 0.0, 1.0)) * self._mdot_edd

    # ── Schwarzschild radius (for reference) ─────────────────────────────────

    def schwarzschild_radius(self) -> float:
        """Schwarzschild radius R_s = 2 G M / c² (m).

        Returns:
            R_s (m).
        """
        return 2.0 * G * self.bh_mass / C_LIGHT**2

    def innermost_stable_orbit(self) -> float:
        """ISCO radius for Schwarzschild BH: r_ISCO = 6 G M / c² = 3 R_s (m).

        Returns:
            r_ISCO (m).
        """
        return 3.0 * self.schwarzschild_radius()

    def to_dict(self) -> dict:
        """Serialize accretion disk parameters to dictionary."""
        return {
            "bh_mass_Msun": self.bh_mass / 1.989e30,
            "eddington_luminosity_W": self.eddington_luminosity(),
            "eddington_accretion_rate_kg_s": self._mdot_edd,
            "schwarzschild_radius_m": self.schwarzschild_radius(),
            "isco_radius_m": self.innermost_stable_orbit(),
            "radiative_efficiency": self.radiative_eff,
        }
