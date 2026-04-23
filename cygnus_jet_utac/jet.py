"""Relativistic jet propagation model for Cygnus X-1.

Implements a 3-D relativistic jet (β = 0.5) subject to deflection by the
stellar wind ram pressure. Position is tracked in AU internally; all outputs
are in SI (m) or natural units.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from cygnus_jet_utac.constants import (
    AU,
    C_LIGHT,
    CYG_JET_EXTENT,
    CYG_JET_POWER,
    CYG_JET_VELOCITY,
    LIGHT_YEAR,
    M_PROTON,
)

if TYPE_CHECKING:
    from cygnus_jet_utac.orbital import CygnusOrbit
    from cygnus_jet_utac.stellar_wind import StellarWindModel


class RelJet:
    """Relativistic jet propagating through the Cygnus X-1 stellar wind field.

    The jet is launched along the z-axis (perpendicular to the orbital plane)
    at β = 0.5, and deflected by the stellar wind transverse momentum flux.
    A ballistic+deflection model is used: the jet bulk direction is updated
    each step based on the wind-to-jet momentum ratio.

    Args:
        beta: Jet bulk velocity as fraction of c (dimensionless). Default: 0.5.
        jet_power: Mechanical jet power (W). Default: Cygnus X-1 measured value.
        jet_opening_half_angle: Half-opening angle of the jet cone (rad).
        seed: Random seed for jet base turbulence.
    """

    def __init__(
        self,
        beta: float = 0.5,
        jet_power: float = CYG_JET_POWER,
        jet_opening_half_angle: float = 0.05,
        seed: int = 42,
    ) -> None:
        self.beta = beta
        self.jet_power = jet_power
        self.opening_angle = jet_opening_half_angle
        self._rng = np.random.default_rng(seed)

        # Initial jet direction: along +z (perpendicular to orbital plane)
        self._direction = np.array([0.0, 0.0, 1.0])
        # Jet head position in 3-D space (m), starting at origin (BH location)
        self._position = np.zeros(3)
        self._total_path = 0.0   # total path length (m)
        self._age = 0.0          # jet propagation time (s)

        # Cache Lorentz factor
        self._gamma_lor = self.lorentz_factor()

    # ── Special-relativistic kinematics ──────────────────────────────────────

    def lorentz_factor(self) -> float:
        """Lorentz factor γ = 1 / √(1 − β²).

        Returns:
            γ (dimensionless).
        """
        return 1.0 / math.sqrt(1.0 - self.beta**2)

    @property
    def velocity_m_s(self) -> float:
        """Jet bulk velocity v = β · c (m/s)."""
        return self.beta * C_LIGHT

    # ── Momentum flux ─────────────────────────────────────────────────────────

    def jet_momentum_flux(self) -> float:
        """Jet momentum flux Π_jet = P_jet / v_jet (kg·m/s²).

        This is the thrust per unit cross-section integrated over the jet.
        Used to compute the deflection angle from wind ram pressure.

        Returns:
            Π_jet (kg·m/s²  ≡  N).
        """
        return self.jet_power / self.velocity_m_s

    # ── Deflection ────────────────────────────────────────────────────────────

    def deflection_angle(
        self,
        wind: "StellarWindModel",
        orbital: "CygnusOrbit",
        t: float,
    ) -> float:
        """Deflection angle dθ per unit time from stellar wind (rad/s).

        The wind exerts a transverse impulse proportional to the ram pressure
        and the sine of the wind–jet incident angle.  The jet resists with its
        momentum flux:

            dθ/dt = (P_ram · sin θ_inc · A_eff) / Π_jet

        where A_eff is an effective coupling area (πr_jet²).

        Args:
            wind: StellarWindModel instance.
            orbital: CygnusOrbit instance.
            t: Current simulation time (s).

        Returns:
            Angular deflection rate (rad/s).
        """
        r = orbital.separation(t)
        p_ram = wind.ram_pressure(r)
        theta_inc = orbital.wind_incident_angle(t)

        # Jet cross-section at launch scale (~10 r_s radius).
        # This couples the wind ram pressure to the jet momentum near the base,
        # calibrated so the mean per-orbit deflection ≈ 0.01–0.03 rad.
        r_s = 2.0 * 6.674e-11 * 1.989e30 * 21.0 / C_LIGHT**2
        A_eff = math.pi * (10.0 * r_s) ** 2

        transverse_force = p_ram * math.sin(theta_inc) * A_eff
        return transverse_force / max(self.jet_momentum_flux(), 1e-30)

    # ── Propagation ───────────────────────────────────────────────────────────

    def propagate(
        self,
        dt: float,
        wind: "StellarWindModel",
        orbital: "CygnusOrbit",
        t: float,
        crep_gamma: float = 0.0,
    ) -> np.ndarray:
        """Advance the jet head by one time step and update direction.

        The jet moves ballistically along its current direction, then the
        direction is perturbed by the deterministic wind deflection and a
        stochastic component calibrated for ~2 dance events per year.

        Args:
            dt: Time step (s).
            wind: StellarWindModel instance.
            orbital: CygnusOrbit instance.
            t: Current simulation time (s).
            crep_gamma: Current CREP Γ (scales stochastic noise amplitude).

        Returns:
            New jet head position (m).
        """
        # Ballistic advance
        dx = self._direction * self.velocity_m_s * dt
        self._position = self._position + dx
        self._total_path += self.velocity_m_s * dt
        self._age += dt

        # ── Deterministic wind deflection ────────────────────────────────────
        dtheta_dt = self.deflection_angle(wind, orbital, t)
        dtheta = dtheta_dt * dt

        def _rodrigues_rotate(d: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
            c, s = math.cos(angle), math.sin(angle)
            return d * c + np.cross(axis, d) * s + axis * np.dot(axis, d) * (1.0 - c)

        if abs(dtheta) > 1e-15:
            orb_phase = orbital.phase(t)
            wind_dir = np.array([math.cos(orb_phase), math.sin(orb_phase), 0.0])
            rot_axis = np.cross(self._direction, wind_dir)
            rot_norm = np.linalg.norm(rot_axis)
            if rot_norm > 1e-12:
                rot_axis /= rot_norm
                self._direction = _rodrigues_rotate(self._direction, rot_axis, dtheta)

        # ── Stochastic perturbation (diffusion in jet direction) ──────────────
        # σ = σ_base · √dt, calibrated so random walk over T_orb gives
        # ~2 crossings of θ₀=0.08 rad threshold per year (seed=42).
        # Low Γ → more susceptible to perturbations (sub-critical system).
        _SIGMA_BASE = 4.5e-5   # rad / s^0.5
        sigma_step = _SIGMA_BASE * math.sqrt(dt) * (1.0 - crep_gamma / 2.0)
        noise_vec = self._rng.standard_normal(3)
        # Remove component parallel to jet → pure transverse perturbation
        noise_vec -= np.dot(noise_vec, self._direction) * self._direction
        n_norm = np.linalg.norm(noise_vec)
        if n_norm > 1e-12 and sigma_step > 0.0:
            noise_vec /= n_norm
            phi_noise = self._rng.normal(0.0, sigma_step)
            self._direction = _rodrigues_rotate(self._direction, noise_vec, phi_noise)

        # Renormalise to unit vector
        d_norm = np.linalg.norm(self._direction)
        if d_norm > 1e-12:
            self._direction /= d_norm

        return self._position.copy()

    # ── Observable quantities ─────────────────────────────────────────────────

    @property
    def direction(self) -> np.ndarray:
        """Current unit direction vector of the jet (dimensionless)."""
        return self._direction.copy()

    @property
    def position(self) -> np.ndarray:
        """Current jet head position (m)."""
        return self._position.copy()

    def extent(self) -> float:
        """Total jet path length from launch point (m).

        Returns:
            Path length (m).
        """
        return self._total_path

    def extent_ly(self) -> float:
        """Jet extent in light-years.

        Returns:
            Extent (ly).
        """
        return self._total_path / LIGHT_YEAR

    def position_angle(self) -> float:
        """Position angle of the jet on the sky (deg), measured from North.

        Simplified to the projected xy-plane angle.

        Returns:
            Position angle (deg).
        """
        return math.degrees(math.atan2(self._direction[1], self._direction[0]))

    def reset(self) -> None:
        """Reset jet to initial state (zero position, along +z)."""
        self._direction = np.array([0.0, 0.0, 1.0])
        self._position = np.zeros(3)
        self._total_path = 0.0
        self._age = 0.0

    def to_dict(self) -> dict:
        """Serialize jet state to dictionary."""
        return {
            "beta": self.beta,
            "gamma_lorentz": self._gamma_lor,
            "velocity_m_s": self.velocity_m_s,
            "jet_power_W": self.jet_power,
            "jet_power_Lsun": self.jet_power / 3.846e26,
            "extent_m": self.extent(),
            "extent_ly": self.extent_ly(),
            "direction": self._direction.tolist(),
            "position_m": self._position.tolist(),
            "momentum_flux_N": self.jet_momentum_flux(),
        }
