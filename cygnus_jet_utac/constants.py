"""Physical constants and Cygnus X-1 system parameters.

All values in SI units unless noted. System parameters from Prabu et al. (2026,
Nature Astronomy, DOI: 10.1038/s41550-026-02828-3).
"""

# ── Fundamental physical constants ────────────────────────────────────────────
C_LIGHT: float = 2.998e8          # m/s  — speed of light
L_SUN: float = 3.846e26           # W    — solar luminosity
M_SUN: float = 1.989e30           # kg   — solar mass
AU: float = 1.496e11              # m    — astronomical unit
LIGHT_YEAR: float = 9.461e15      # m    — light-year
G: float = 6.674e-11              # m³/(kg·s²)  — gravitational constant
SIGMA_T: float = 6.652e-29        # m²   — Thomson cross-section
M_PROTON: float = 1.673e-27       # kg   — proton mass
K_BOLTZMANN: float = 1.381e-23    # J/K  — Boltzmann constant

# ── Cygnus X-1 system parameters (Prabu et al. 2026) ─────────────────────────
CYG_BH_MASS: float = 21.0 * M_SUN           # kg  — black hole mass
CYG_COMPANION_MASS: float = 41.0 * M_SUN    # kg  — blue supergiant (HDE 226868)
CYG_ORBITAL_PERIOD: float = 5.6 * 86400.0   # s   — orbital period
CYG_DISTANCE: float = 7000.0 * LIGHT_YEAR   # m   — for angular scaling only
CYG_JET_POWER: float = 10_000.0 * L_SUN     # W   — measured jet power
CYG_JET_VELOCITY: float = 0.5 * C_LIGHT     # m/s — measured jet bulk velocity
CYG_JET_EXTENT: float = 16.0 * LIGHT_YEAR   # m   — observed radio jet extent
CYG_ACCRETION_EFFICIENCY: float = 0.10      # dimensionless — η = P_jet/P_acc
CYG_OBSERVATION_YEARS: float = 18.0         # yr  — VLBI baseline (2006–2024)

# ── Stellar wind parameters (HDE 226868 OB supergiant) ────────────────────────
# Terminal wind velocity (CAK model): v_inf ≈ 2100 km/s for OB supergiants.
# In the CREP framework this is scaled to 3×CYG_JET_VELOCITY to set the
# CREP Resonance field amplitude (model parameter, not the physical v_inf).
CYG_WIND_VINF: float = 3.0 * CYG_JET_VELOCITY   # m/s — CREP model amplitude
CYG_WIND_BETA: float = 0.8                        # CAK β-velocity law exponent
CYG_WIND_MDOT: float = 2.5e-6 * (M_SUN / (365.25 * 86400))  # kg/s ≈ 2.5e-6 M☉/yr

# Stellar radius from mass-luminosity relation for OB supergiants (R ~ 20–25 R☉)
R_SUN: float = 6.957e8                        # m
CYG_COMPANION_RADIUS: float = 22.0 * R_SUN   # m

# ── Accretion disk parameters ─────────────────────────────────────────────────
ACCRETION_RADIATIVE_EFF: float = 0.10    # η_rad — standard thin disk efficiency

# ── UTAC calibrated defaults (GenesisAeon ERA5 baseline) ─────────────────────
UTAC_R_DEFAULT: float = 0.12     # intrinsic UTAC growth rate
UTAC_SIGMA_DEFAULT: float = 2.2  # CREP–UTAC coupling coefficient σ
SIGMA_PHI: float = 1.0 / 16.0   # Frame Principle minimum variance ratio

# ── Derived quantities ────────────────────────────────────────────────────────
CYG_ORBITAL_ANGULAR_VEL: float = 2.0 * 3.14159265358979 / CYG_ORBITAL_PERIOD  # rad/s

# Keplerian semi-major axis a from Kepler's third law: a³ = G(M1+M2)T²/(4π²)
import math as _math  # noqa: E402

CYG_ORBITAL_SEPARATION: float = (
    G * (CYG_BH_MASS + CYG_COMPANION_MASS) * CYG_ORBITAL_PERIOD**2
    / (4.0 * _math.pi**2)
) ** (1.0 / 3.0)   # m

del _math
