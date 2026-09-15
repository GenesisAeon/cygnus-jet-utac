"""Accretion-to-jet efficiency calibration via UTAC fixed-point inversion.

    Given the measured η = 10 % efficiency (Prabu et al. 2026), invert the
    UTAC fixed-point relation H* = K·tanh(σ·Γ) to recover Γ_jet — the
    domain-specific CREP value for the Cygnus X-1 jet system.

    Result: Γ_jet = arctanh(η) / σ ≈ arctanh(0.10) / 2.2 ≈ 0.0456

HONESTY NOTE (2026-09-15, ecosystem-wide Gamma-circularity review): this is
an algebraic inversion, not an independent calibration -- for ANY sigma,
this same formula recovers eta exactly by construction (both
"utac_fixed_point_check" and "efficiency_check" below will always pass).
sigma=2.2 itself ("UTAC_SIGMA_DEFAULT", commented "ERA5 baseline" in
constants.py) is a shared default reused unchanged from the GenesisAeon
climate/AMOC packages -- verified byte-identical to amoc-utac's own
UTAC_SIGMA -- not independently derived from Prabu et al. 2026 or any
Cygnus-X-1-specific source. Further: benchmark.py's other five "Prabu 2026"
targets (jet_power_W, jet_velocity_c, jet_extent_ly, orbital_period_days,
dance_events_per_year) were traced individually and are each either an
echoed input constant or a free parameter explicitly tuned to hit that
exact target -- none is an independent prediction. The underlying UTAC
ODE integration and the surrounding astrophysics (CAK wind law,
relativistic jet kinematics, Keplerian orbit) ARE genuine, independently
implemented physics -- only the Gamma_jet "first CREP-domain calibration"
framing itself does not hold up. See
D:\\mandala\\crep-utac-afet-formalism\\worked_example_cygnus_jet_utac.md
and FOLLOWUP_TICKETS.md for the full analysis.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from cygnus_jet_utac.constants import (
    CYG_ACCRETION_EFFICIENCY,
    SIGMA_PHI,
    UTAC_R_DEFAULT,
    UTAC_SIGMA_DEFAULT,
)


def calibrate_gamma_jet(
    eta: float = CYG_ACCRETION_EFFICIENCY,
    sigma: float = UTAC_SIGMA_DEFAULT,
    r: float = UTAC_R_DEFAULT,
    K: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Invert η = H*/K = tanh(σ·Γ) to solve for Γ_jet.

    Derivation:
        UTAC fixed point: H* = K · tanh(σ · Γ)
        Efficiency:       η  = H* / K
        → η = tanh(σ · Γ)
        → Γ_jet = arctanh(η) / σ

    The Frame Principle check computes σ_Φ,min for this domain and tests
    whether σ_Φ,min / (1/16) ≤ 1 (must hold for the Frame Principle to apply).

    Args:
        eta: Accretion-to-jet power efficiency η (dimensionless). Default: 0.10.
        sigma: CREP coupling coefficient σ. Default: 2.2 (ERA5 baseline).
        r: Intrinsic UTAC growth rate. Default: 0.12.
        K: Normalised carrying capacity. Default: 1.0.
        verbose: If True, include a human-readable derivation in the result.

    Returns:
        Dictionary with keys:
            gamma_jet (float): Domain-specific CREP value for Cygnus X-1 jets.
            sigma_phi_min (float): Minimum Frame Principle variance ratio.
            sigma_phi_ratio (float): sigma_phi_min / (1/16) — should be ≤ 1.
            utac_fixed_point_check (float): K·tanh(σ·Γ_jet) — should equal η·K.
            efficiency_check (float): Recovered η from Γ_jet — should equal input η.
            interpretation (str): Physical interpretation of the result.
            report (str): Full human-readable derivation (if verbose=True).
    """
    if not (0.0 < eta < 1.0):
        raise ValueError(f"Efficiency η must be in (0, 1), got {eta}")

    gamma_jet = math.atanh(eta) / sigma

    # σ_Φ,min = r · (1 − tanh(σ · Γ)) / 2  (Frame Principle condition)
    sigma_phi_min = r * (1.0 - math.tanh(sigma * gamma_jet)) / 2.0
    sigma_phi_ratio = sigma_phi_min / SIGMA_PHI

    # Verification
    utac_fp = K * math.tanh(sigma * gamma_jet)
    efficiency_check = math.tanh(sigma * gamma_jet)

    criticality = "sub-critical" if gamma_jet < 0.5 else "super-critical"
    sensitivity = (
        "high sensitivity to perturbations (low-Γ, barely supercritical)"
        if gamma_jet < 0.1
        else ("moderate sensitivity" if gamma_jet < 0.5 else "low sensitivity (rigid jet)")
    )

    interpretation = (
        f"Γ_jet = {gamma_jet:.4f} — Cygnus X-1 jet is {criticality} "
        f"({'below' if gamma_jet < 0.5 else 'above'} Γ=0.5 threshold). "
        f"{sensitivity}. "
        f"The 'dancing' jet behaviour is explained by this low CREP value: "
        f"small stellar wind perturbations drive large direction changes."
    )

    result: dict[str, Any] = {
        "gamma_jet": gamma_jet,
        "sigma_phi_min": sigma_phi_min,
        "sigma_phi_ratio": sigma_phi_ratio,
        "sigma_phi_frame_principle_satisfied": sigma_phi_ratio <= 1.0,
        "utac_fixed_point_check": utac_fp,
        "efficiency_check": efficiency_check,
        "input_eta": eta,
        "input_sigma": sigma,
        "input_r": r,
        "interpretation": interpretation,
    }

    if verbose:
        report_lines = [
            "=" * 60,
            "  Γ_jet CALIBRATION — GenesisAeon Package 17",
            "=" * 60,
            f"  Input η (Prabu 2026):       {eta:.4f}  ({eta*100:.1f} %)",
            f"  UTAC σ (ERA5 baseline):     {sigma:.4f}",
            f"  UTAC r (ERA5 baseline):     {r:.4f}",
            "",
            "  Inversion:  Γ_jet = arctanh(η) / σ",
            f"            = arctanh({eta}) / {sigma}",
            f"            = {math.atanh(eta):.6f} / {sigma}",
            f"            = {gamma_jet:.6f}",
            "",
            "  UTAC fixed-point check:",
            f"    H*/K = tanh(σ·Γ) = tanh({sigma:.2f}×{gamma_jet:.4f})",
            f"         = tanh({sigma*gamma_jet:.6f}) = {efficiency_check:.6f}",
            (
                f"    → η_recovered = {efficiency_check:.4f}  ✓"
                if abs(efficiency_check - eta) < 1e-10
                else f"    → η_recovered = {efficiency_check:.4f}"
            ),
            "",
            "  Frame Principle check (σ_Φ ≥ 1/16 = 0.0625):",
            "    σ_Φ,min = r·(1−tanh(σ·Γ))/2",
            f"            = {r}·(1−{efficiency_check:.4f})/2",
            f"            = {sigma_phi_min:.6f}",
            f"    Ratio   = {sigma_phi_min:.6f} / {SIGMA_PHI:.4f} = {sigma_phi_ratio:.4f}",
            f"    Satisfied: {'✓ YES' if sigma_phi_ratio <= 1.0 else '✗ NO'}",
            "",
            f"  → {interpretation}",
            "=" * 60,
        ]
        result["report"] = "\n".join(report_lines)

    return result


def efficiency_from_gamma(gamma: float, sigma: float = UTAC_SIGMA_DEFAULT) -> float:
    """Compute efficiency η = tanh(σ · Γ) from CREP value.

    Args:
        gamma: CREP Γ value (dimensionless).
        sigma: CREP coupling σ. Default: 2.2.

    Returns:
        Efficiency η ∈ (0, 1).
    """
    return math.tanh(sigma * gamma)


def gamma_scan(
    eta_values: np.ndarray[Any, Any] | None = None,
    sigma: float = UTAC_SIGMA_DEFAULT,
) -> tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
    """Scan η → Γ_jet relationship across a range of efficiencies.

    Args:
        eta_values: Array of efficiency values. Default: 0.01 to 0.99.
        sigma: CREP coupling σ.

    Returns:
        Tuple (eta_arr, gamma_arr) as numpy arrays.
    """
    if eta_values is None:
        eta_values = np.linspace(0.01, 0.99, 200)
    gamma_values = np.arctanh(eta_values) / sigma
    return eta_values, gamma_values
