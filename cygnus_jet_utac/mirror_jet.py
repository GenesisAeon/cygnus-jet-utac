"""Mirror-Machine for detecting phase transitions in jet orientation.

Applies the GenesisAeon Mirror-Machine concept: the current jet direction
d(t) is compared against the τ-delayed mirror state d(t − τ), where
τ = one orbital period. A phase transition ("jet dance event") is triggered
when the Jensen-Shannon divergence D_mirror exceeds the CREP-adaptive
threshold θ_PT = θ₀ · (1 − Γ(t)/2).
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np

from cygnus_jet_utac.constants import CYG_ORBITAL_PERIOD


@dataclass
class DanceEvent:
    """A detected jet direction change (phase transition) event.

    Args:
        t: Simulation time of the event (s).
        t_years: Simulation time in years.
        divergence: D_mirror at the moment of detection.
        threshold: θ_PT threshold at the moment of detection.
        direction_before: Jet direction at t − τ.
        direction_after: Jet direction at t.
        deflection_angle_deg: Angular change between the two directions (deg).
        crep_gamma: CREP Γ at the time of the event.
    """
    t: float
    t_years: float
    divergence: float
    threshold: float
    direction_before: list
    direction_after: list
    deflection_angle_deg: float
    crep_gamma: float
    metadata: dict = field(default_factory=dict)


class JetMirrorMachine:
    """Mirror-Machine for detecting phase transitions in jet orientation.

    The jet direction history is stored in a ring buffer of length N_buf.
    The buffer is indexed so that d(t − τ) is retrieved by stepping back
    τ / dt samples in the buffer.

    Phase events are triggered when:

        D_mirror(t) = angular_divergence(d(t), d(t−τ)) > θ_PT(t)
        θ_PT(t) = θ₀ · (1 − Γ(t) / 2)

    A lower Γ (sub-critical system, like Cygnus X-1 at Γ ≈ 0.046) gives a
    HIGHER threshold, meaning events are less likely — consistent with the
    jet being highly responsive to external perturbations but not constantly
    triggering.  Physically, a low-Γ jet "dances" when perturbed significantly.

    Args:
        tau: Mirror delay — one orbital period (s).
        dt: Simulation time step (s).
        theta0: Base phase-transition threshold (rad).
        gamma0: Reference CREP value for threshold scaling.
    """

    def __init__(
        self,
        tau: float = CYG_ORBITAL_PERIOD,
        dt: float = 3600.0,
        theta0: float = 0.08,
        gamma0: float = 0.5,
    ) -> None:
        self.tau = tau
        self.dt = dt
        self.theta0 = theta0
        self.gamma0 = gamma0

        # Number of steps in one mirror delay
        self._n_delay = max(1, int(round(tau / dt)))
        # Ring buffer for direction history
        self._buf: deque[np.ndarray] = deque(
            [np.array([0.0, 0.0, 1.0])] * (self._n_delay + 1),
            maxlen=self._n_delay + 1,
        )
        self._dance_events: list[DanceEvent] = []
        self._last_event_t: float = -tau  # prevent back-to-back triggers
        self._current_divergence: float = 0.0
        self._current_threshold: float = theta0
        self._n_steps: int = 0  # step counter for warmup guard

    # ── Update ────────────────────────────────────────────────────────────────

    def update(
        self,
        direction: np.ndarray,
        crep_gamma: float,
        t: float,
    ) -> bool:
        """Update the Mirror-Machine with a new jet direction sample.

        Args:
            direction: Current jet unit direction vector (3-D).
            crep_gamma: Current aggregate CREP Γ value (dimensionless).
            t: Current simulation time (s).

        Returns:
            True if a phase transition event was triggered, False otherwise.
        """
        d_now = np.asarray(direction, dtype=float)
        norm = np.linalg.norm(d_now)
        if norm > 1e-12:
            d_now = d_now / norm

        # Retrieve the τ-delayed direction (oldest in the buffer)
        d_delayed = self._buf[0]

        # Jensen-Shannon-like angular divergence: D = arccos(|d·d_τ|)
        dot = float(np.clip(np.dot(d_now, d_delayed), -1.0, 1.0))
        self._current_divergence = math.acos(abs(dot))  # rad

        # Adaptive threshold: higher Gamma → lower threshold (subtler transitions counted).
        self._current_threshold = self.theta0 * (1.0 - crep_gamma / 2.0)

        # Push current direction into buffer
        self._buf.append(d_now.copy())
        self._n_steps += 1

        # Only trigger after the buffer is fully populated (warmup = n_delay steps)
        # and the refractory period (0.5 τ) has elapsed since the last event.
        triggered = (
            self._n_steps > self._n_delay
            and self._current_divergence > self._current_threshold
            and (t - self._last_event_t) > 0.5 * self.tau
        )

        if triggered:
            deflection_deg = math.degrees(self._current_divergence)
            event = DanceEvent(
                t=t,
                t_years=t / (365.25 * 86400.0),
                divergence=self._current_divergence,
                threshold=self._current_threshold,
                direction_before=d_delayed.tolist(),
                direction_after=d_now.tolist(),
                deflection_angle_deg=deflection_deg,
                crep_gamma=crep_gamma,
                metadata={"refractory_delay_tau": (t - self._last_event_t) / self.tau},
            )
            self._dance_events.append(event)
            self._last_event_t = t

        return triggered

    # ── Accessors ─────────────────────────────────────────────────────────────

    def divergence(self) -> float:
        """Current angular divergence D_mirror (rad)."""
        return self._current_divergence

    def threshold(self) -> float:
        """Current adaptive threshold θ_PT (rad)."""
        return self._current_threshold

    def get_dance_events(self) -> list[DanceEvent]:
        """Return list of all detected jet dance events."""
        return list(self._dance_events)

    def dance_events_per_year(self, total_years: float) -> float:
        """Average dance event rate (events/year).

        Args:
            total_years: Total simulation duration (years).

        Returns:
            Events per year.
        """
        if total_years <= 0.0:
            return 0.0
        return len(self._dance_events) / total_years

    def reset(self) -> None:
        """Reset internal state (clears events and buffer)."""
        self._buf = deque(
            [np.array([0.0, 0.0, 1.0])] * (self._n_delay + 1),
            maxlen=self._n_delay + 1,
        )
        self._dance_events.clear()
        self._last_event_t = -self.tau
        self._current_divergence = 0.0
        self._n_steps = 0

    def to_dict(self) -> dict:
        """Serialize Mirror-Machine state to dictionary."""
        return {
            "tau_s": self.tau,
            "tau_days": self.tau / 86400.0,
            "dt_s": self.dt,
            "n_delay_steps": self._n_delay,
            "theta0_rad": self.theta0,
            "n_dance_events": len(self._dance_events),
            "current_divergence_rad": self._current_divergence,
            "current_threshold_rad": self._current_threshold,
        }
