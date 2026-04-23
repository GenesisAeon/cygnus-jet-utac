"""Binary orbital dynamics for the Cygnus X-1 system.

Implements a Keplerian circular orbit for the BH + Blue Supergiant binary.
The orbital phase φ(t) drives the time-varying stellar wind incident angle on
the jet, which feeds into the CREP Resonance component.
"""

from __future__ import annotations

import math

import numpy as np

from cygnus_jet_utac.constants import (
    CYG_BH_MASS,
    CYG_COMPANION_MASS,
    CYG_ORBITAL_PERIOD,
    G,
)


class CygnusOrbit:
    """Keplerian circular orbit for the Cygnus X-1 BH + Blue Supergiant binary.

    Cygnus X-1 has a very low eccentricity (e ≈ 0.02), so the orbit is
    approximated as circular. The orbital phase drives the stellar wind
    incident angle on the jet axis, providing the time-varying external
    CREP modulation field.

    Args:
        period: Orbital period (s). Default: 5.6 days.
        bh_mass: Black hole mass (kg).
        companion_mass: Companion (blue supergiant) mass (kg).
        t0: Reference epoch (s) for zero phase.
    """

    def __init__(
        self,
        period: float = CYG_ORBITAL_PERIOD,
        bh_mass: float = CYG_BH_MASS,
        companion_mass: float = CYG_COMPANION_MASS,
        t0: float = 0.0,
    ) -> None:
        self.period = period
        self.bh_mass = bh_mass
        self.companion_mass = companion_mass
        self.t0 = t0
        self._total_mass = bh_mass + companion_mass
        self._omega = 2.0 * math.pi / period  # rad/s
        # Derive separation from Kepler's third law using the provided masses/period
        self._separation = (
            G * self._total_mass * period**2 / (4.0 * math.pi**2)
        ) ** (1.0 / 3.0)

    @property
    def separation_m(self) -> float:
        """Mean orbital separation (m) — constant for circular orbit."""
        return self._separation

    @property
    def angular_velocity(self) -> float:
        """Orbital angular velocity ω = 2π/T (rad/s)."""
        return self._omega

    def phase(self, t: float) -> float:
        """Orbital phase φ(t) = ω·(t − t₀) mod 2π.

        Args:
            t: Simulation time (s).

        Returns:
            Phase φ ∈ [0, 2π] (rad).
        """
        return (self._omega * (t - self.t0)) % (2.0 * math.pi)

    def separation(self, t: float) -> float:
        """Orbital separation between BH and companion (m).

        Constant for circular orbit; provided for API completeness.

        Args:
            t: Simulation time (s).

        Returns:
            Separation (m).
        """
        return self._separation

    def position_bh(self, t: float) -> np.ndarray:
        """3-D position of the BH in the centre-of-mass frame (m).

        Args:
            t: Simulation time (s).

        Returns:
            Position vector [x, y, 0] (m), orbit in xy-plane.
        """
        phi = self.phase(t)
        # BH orbits at distance r_bh = a * M2/(M1+M2) from CoM
        r_bh = self._separation * self.companion_mass / self._total_mass
        return np.array([r_bh * math.cos(phi), r_bh * math.sin(phi), 0.0])

    def position_companion(self, t: float) -> np.ndarray:
        """3-D position of the Blue Supergiant in the CoM frame (m).

        Args:
            t: Simulation time (s).

        Returns:
            Position vector [x, y, 0] (m).
        """
        phi = self.phase(t)
        r_comp = self._separation * self.bh_mass / self._total_mass
        return np.array([-r_comp * math.cos(phi), -r_comp * math.sin(phi), 0.0])

    def wind_incident_angle(self, t: float) -> float:
        """Angle between the stellar wind direction and the jet axis (rad).

        The jet is launched perpendicular to the orbital plane (z-direction).
        The stellar wind blows radially outward from the companion, so the
        incident angle on the jet is the angle between the companion-to-BH
        vector and the jet axis (ẑ).

        For a circular orbit the incident angle oscillates between 0 and π/2,
        modulated by orbital phase, creating the CREP Resonance oscillation.

        Args:
            t: Simulation time (s).

        Returns:
            Incident angle θ ∈ [0, π/2] (rad).
        """
        phi = self.phase(t)
        # Wind-jet coupling is strongest at conjunction (φ = π/2, 3π/2)
        # and weakest at quadrature (φ = 0, π).
        # The effective projected coupling is |sin(φ)| modulated by geometry.
        return math.pi / 4.0 * (1.0 + math.sin(phi))

    def orbital_velocity(self) -> float:
        """Circular orbital velocity of the BH (m/s).

        Returns:
            v_orb = ω · r_bh (m/s).
        """
        r_bh = self._separation * self.companion_mass / self._total_mass
        return self._omega * r_bh

    def roche_lobe_radius(self) -> float:
        """Approximate Roche lobe radius of the BH (Eggleton 1983, m).

        Returns:
            R_L (m).
        """
        q = self.bh_mass / self.companion_mass
        numer = 0.49 * q ** (2.0 / 3.0)
        denom = 0.6 * q ** (2.0 / 3.0) + math.log(1.0 + q ** (1.0 / 3.0))
        return self._separation * numer / denom

    def to_dict(self) -> dict:
        """Serialize orbital parameters to dictionary."""
        return {
            "period_s": self.period,
            "period_days": self.period / 86400.0,
            "separation_AU": self._separation / 1.496e11,
            "angular_velocity_rad_s": self._omega,
            "orbital_velocity_m_s": self.orbital_velocity(),
            "bh_mass_Msun": self.bh_mass / 1.989e30,
            "companion_mass_Msun": self.companion_mass / 1.989e30,
        }
