"""
Lightweight stubs for genesis-os core imports.

# STUB — replace with genesis.core when genesis-os is installed:
#   from genesis.core.utac import UTAC_ODE, UTACParams
#   from genesis.core.crep import CREPTensor
#   from genesis.mirror.phase_loop import PhaseTransitionLoop
#   from genesis.core.lagrangian import UnifiedLagrangian
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# ── UTAC ──────────────────────────────────────────────────────────────────────

@dataclass
class UTACParams:
    """Parameters for the UTAC ODE system.

    Args:
        r: Intrinsic growth rate (dimensionless).
        sigma: CREP–UTAC coupling coefficient σ (dimensionless).
        K: Carrying capacity (normalised units, default 1.0).
        noise_scale: Stochastic noise amplitude (dimensionless).
    """
    r: float = 0.12
    sigma: float = 2.2
    K: float = 1.0
    noise_scale: float = 0.0


class UTAC_ODE:
    """Unified Theory of Adaptive Criticality — ODE integrator.

    Implements the UTAC logistic equation with CREP-adaptive carrying capacity:

        dH/dt = r · H · (1 − H / H*)
        H*(t) = K · tanh(σ · Γ(t))

    The system converges to H* as t → ∞ for any Γ > 0.
    """

    def __init__(self, params: UTACParams, H0: float = 0.01) -> None:
        self.params = params
        self._H: float = H0
        self._t: float = 0.0

    @property
    def H(self) -> float:
        """Current state variable H (normalised accretion rate)."""
        return self._H

    @property
    def t(self) -> float:
        """Current simulation time (s)."""
        return self._t

    def fixed_point(self, gamma: float) -> float:
        """Compute UTAC fixed point H* = K·tanh(σ·Γ).

        Args:
            gamma: Current CREP Γ value (dimensionless).

        Returns:
            H* in normalised units [0, K].
        """
        return self.params.K * math.tanh(self.params.sigma * gamma)

    def derivative(self, H: float, gamma: float) -> float:
        """Compute dH/dt for given state and CREP value.

        Args:
            H: Current state variable (normalised).
            gamma: Current CREP Γ (dimensionless).

        Returns:
            dH/dt (s⁻¹).
        """
        H_star = self.fixed_point(gamma)
        if H_star <= 0.0:
            return 0.0
        return self.params.r * H * (1.0 - H / H_star)

    def step(
        self,
        gamma: float,
        dt: float,
        rng: np.random.Generator | None = None,
    ) -> tuple[float, float]:
        """Advance the UTAC ODE by one time step using 4th-order Runge-Kutta.

        Args:
            gamma: CREP Γ at current time (dimensionless).
            dt: Time step (s).
            rng: Optional random generator for stochastic noise.

        Returns:
            Tuple of (H_new, dH_dt) — new state and instantaneous derivative.
        """
        p = self.params
        k1 = self.derivative(self._H, gamma)
        k2 = self.derivative(self._H + 0.5 * dt * k1, gamma)
        k3 = self.derivative(self._H + 0.5 * dt * k2, gamma)
        k4 = self.derivative(self._H + dt * k3, gamma)

        dH_dt = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        noise = 0.0
        if rng is not None and p.noise_scale > 0.0:
            noise = p.noise_scale * rng.standard_normal() * math.sqrt(dt)

        self._H = max(0.0, min(p.K, self._H + dt * dH_dt + noise))
        self._t += dt
        return self._H, dH_dt


# ── CREP ──────────────────────────────────────────────────────────────────────

@dataclass
class CREPTensor:
    """CREP (Coherence–Resonance–Emergence–Pattern) tensor.

    Each component ∈ [0, 1]. The aggregate CREP value Γ is the geometric mean:
        Γ = (C · R · E · P)^(1/4)

    Args:
        C: Coherence — temporal phase alignment of accretion and wind.
        R: Resonance — wind ram pressure at resonance frequency.
        E: Emergence — fraction of accretion power channelled to jet.
        P: Pattern — repetition / orbital periodicity strength.
    """
    C: float = 0.5
    R: float = 0.5
    E: float = 0.5
    P: float = 0.5

    @property
    def Gamma(self) -> float:
        """Γ = (C · R · E · P)^(1/4) — aggregate CREP value."""
        return float(
            (max(0.0, self.C) * max(0.0, self.R)
             * max(0.0, self.E) * max(0.0, self.P)) ** 0.25
        )

    def to_dict(self) -> dict[str, float]:
        """Serialize to dictionary for logging and Zenodo export."""
        return {"C": self.C, "R": self.R, "E": self.E, "P": self.P, "Gamma": self.Gamma}

    def clamp(self) -> CREPTensor:
        """Return a copy with all components clamped to [0, 1]."""
        return CREPTensor(
            C=max(0.0, min(1.0, self.C)),
            R=max(0.0, min(1.0, self.R)),
            E=max(0.0, min(1.0, self.E)),
            P=max(0.0, min(1.0, self.P)),
        )


# ── Phase Transition ──────────────────────────────────────────────────────────

@dataclass
class PhaseTransitionEvent:
    """A single recorded phase transition event.

    Args:
        t: Simulation time of the event (s).
        state: Snapshot of UTAC/CREP state at transition.
        metadata: Domain-specific labels (e.g. jet deflection angle).
    """
    t: float
    state: dict
    metadata: dict


class PhaseTransitionLoop:
    """Phase transition event logger and dispatcher."""

    def __init__(self) -> None:
        self._events: list[PhaseTransitionEvent] = []

    def record_event(
        self, t: float, state: dict, metadata: dict
    ) -> PhaseTransitionEvent:
        """Record a phase transition event.

        Args:
            t: Simulation time (s).
            state: UTAC/CREP state snapshot.
            metadata: Additional domain-specific data.

        Returns:
            The recorded PhaseTransitionEvent.
        """
        event = PhaseTransitionEvent(t=t, state=state, metadata=metadata)
        self._events.append(event)
        return event

    def get_events(self) -> list[PhaseTransitionEvent]:
        """Return all recorded events (copy)."""
        return list(self._events)

    def clear(self) -> None:
        """Clear all recorded events."""
        self._events.clear()

    def __len__(self) -> int:
        return len(self._events)


# ── Lagrangian ────────────────────────────────────────────────────────────────

class UnifiedLagrangian:
    """Unified Lagrangian — Newtonian gravitational potential V(r).

    Args:
        M1: Mass of body 1 (kg).
        M2: Mass of body 2 (kg).
        G_const: Gravitational constant (m³·kg⁻¹·s⁻²).
    """

    def __init__(self, M1: float, M2: float, G_const: float = 6.674e-11) -> None:
        self.M1 = M1
        self.M2 = M2
        self.G = G_const
        self.mu = M1 * M2 / (M1 + M2)

    def potential(self, r: float) -> float:
        """Gravitational potential V(r) = −G·M₁·M₂/r.

        Args:
            r: Separation between bodies (m).

        Returns:
            V(r) in joules (J).
        """
        return -self.G * self.M1 * self.M2 / r

    def gradient(self, r: float) -> float:
        """Radial force magnitude |F(r)| = G·M₁·M₂/r².

        Args:
            r: Separation (m).

        Returns:
            Force magnitude (N).
        """
        return self.G * self.M1 * self.M2 / r**2
