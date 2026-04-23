"""Blue Supergiant stellar wind model — external CREP field modulator.

Implements a CAK (Castor, Abbott & Klein 1975) line-driven wind with a
β-velocity law. The wind ram pressure at the jet location becomes the
Resonance (R) component of the CREP tensor, coupling the orbital dynamics
to the UTAC state of the jet.
"""

from __future__ import annotations

import math

import numpy as np

from cygnus_jet_utac.constants import (
    CYG_COMPANION_MASS,
    CYG_COMPANION_RADIUS,
    CYG_ORBITAL_SEPARATION,
    CYG_WIND_BETA,
    CYG_WIND_MDOT,
    CYG_WIND_VINF,
    CYG_ORBITAL_PERIOD,
)


class StellarWindModel:
    """CAK β-velocity law stellar wind for the HDE 226868 Blue Supergiant.

    The wind terminal velocity parameter ``v_inf`` is set to
    ``3 × CYG_JET_VELOCITY`` as the CREP Resonance amplitude calibration
    (see constants.py). This gives the correct ram-pressure ratio for the
    CREP R-component.

    The wind ram pressure at the jet location is mapped to R ∈ [0, 1] and
    is maximised once per orbital period — the coherence resonance condition.

    Args:
        v_inf: Terminal wind velocity (m/s). Default: 3 × v_jet.
        beta: β-exponent for velocity law (dimensionless). Default: 0.8.
        mdot: Mass-loss rate (kg/s).
        R_star: Stellar radius (m).
        orbital_period: Used to compute resonance modulation (s).
    """

    def __init__(
        self,
        v_inf: float = CYG_WIND_VINF,
        beta: float = CYG_WIND_BETA,
        mdot: float = CYG_WIND_MDOT,
        R_star: float = CYG_COMPANION_RADIUS,
        orbital_period: float = CYG_ORBITAL_PERIOD,
    ) -> None:
        self.v_inf = v_inf
        self.beta = beta
        self.mdot = mdot
        self.R_star = R_star
        self.orbital_period = orbital_period
        # Reference pressure at the mean orbital separation
        self._p_ref = self._ram_pressure_raw(CYG_ORBITAL_SEPARATION)

    # ── Velocity profile ──────────────────────────────────────────────────────

    def velocity(self, r: float) -> float:
        """CAK β-velocity law: v(r) = v_inf · (1 − R*/r)^β.

        Args:
            r: Distance from stellar centre (m).  Must satisfy r ≥ R_star.

        Returns:
            Wind velocity (m/s).
        """
        if r <= self.R_star:
            return 0.0
        return self.v_inf * (1.0 - self.R_star / r) ** self.beta

    # ── Density ───────────────────────────────────────────────────────────────

    def number_density(self, r: float) -> float:
        """Wind proton number density from mass-continuity (m⁻³).

        n(r) = Ṁ / (4π r² m_p v(r))

        Args:
            r: Distance from stellar centre (m).

        Returns:
            Number density (m⁻³).
        """
        from cygnus_jet_utac.constants import M_PROTON
        v = self.velocity(r)
        if v <= 0.0:
            return 0.0
        return self.mdot / (4.0 * math.pi * r**2 * M_PROTON * v)

    def mass_density(self, r: float) -> float:
        """Wind mass density (kg/m³).

        ρ(r) = Ṁ / (4π r² v(r))

        Args:
            r: Distance from stellar centre (m).

        Returns:
            Mass density (kg/m³).
        """
        v = self.velocity(r)
        if v <= 0.0:
            return 0.0
        return self.mdot / (4.0 * math.pi * r**2 * v)

    # ── Ram pressure ──────────────────────────────────────────────────────────

    def _ram_pressure_raw(self, r: float) -> float:
        """Wind dynamic (ram) pressure P_ram = ½ ρ v² (Pa).

        Args:
            r: Distance from stellar centre (m).

        Returns:
            Ram pressure (Pa).
        """
        rho = self.mass_density(r)
        v = self.velocity(r)
        return 0.5 * rho * v**2

    def ram_pressure(self, r: float) -> float:
        """Wind ram pressure at distance r (Pa).

        Args:
            r: Distance from stellar centre (m).

        Returns:
            Ram pressure (Pa).
        """
        return self._ram_pressure_raw(r)

    # ── CREP R-component ──────────────────────────────────────────────────────

    def crep_R_component(self, t: float, jet_pos: np.ndarray) -> float:
        """Map wind ram pressure to CREP Resonance component R ∈ [0, 1].

        R is maximised when the wind–jet coupling is at the resonance
        frequency (once per orbital period — the coherence resonance
        condition).  The resonance modulation is:

            R(t) = R_base · (1 + cos(2πt/T_orb)) / 2

        where R_base = P_ram(r) / P_ram(a_orb) is the normalised pressure.

        Args:
            t: Simulation time (s).
            jet_pos: Jet axis reference position [x, y, z] (m).

        Returns:
            R ∈ [0, 1] (dimensionless).
        """
        # The wind interacts with the jet at its LAUNCH SCALE (near the BH),
        # not at the jet head which is far away. Use the orbital separation as
        # the relevant interaction radius — the wind-jet coupling scale.
        r = CYG_ORBITAL_SEPARATION

        p_ram = self._ram_pressure_raw(r)
        if self._p_ref <= 0.0:
            return 0.0

        # Normalised base resonance level
        R_base = min(1.0, p_ram / self._p_ref)

        # Orbital resonance modulation — peaks once per orbital period
        resonance_phase = 2.0 * math.pi * t / self.orbital_period
        modulation = 0.5 * (1.0 + math.cos(resonance_phase))

        return float(np.clip(R_base * (0.3 + 0.7 * modulation), 0.0, 1.0))

    # ── Wind force on jet ─────────────────────────────────────────────────────

    def wind_force_on_jet(
        self, r: float, jet_direction: np.ndarray, incident_angle: float
    ) -> np.ndarray:
        """Transverse force per unit length exerted by wind on jet (N/m).

        Uses the momentum-flux approach: F/L = P_ram · sin(θ) · ĵ_perp,
        where θ is the wind–jet incident angle and ĵ_perp is the unit vector
        perpendicular to the jet lying in the wind–jet plane.

        Args:
            r: Separation from companion to jet (m).
            jet_direction: Unit vector along jet axis (dimensionless).
            incident_angle: Wind–jet angle θ (rad).

        Returns:
            Transverse force per unit length (N/m).
        """
        p_ram = self.ram_pressure(r)
        # Perpendicular to jet in the xy-plane
        perp = np.array([-jet_direction[1], jet_direction[0], 0.0])
        norm = np.linalg.norm(perp)
        if norm < 1e-12:
            return np.zeros(3)
        perp /= norm
        return p_ram * math.sin(incident_angle) * perp

    def to_dict(self) -> dict:
        """Serialize wind parameters to dictionary."""
        return {
            "v_inf_m_s": self.v_inf,
            "beta": self.beta,
            "mdot_Msun_yr": self.mdot * 365.25 * 86400 / 1.989e30,
            "R_star_Rsun": self.R_star / 6.957e8,
            "orbital_period_days": self.orbital_period / 86400.0,
        }
