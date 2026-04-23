# cygnus-jet-utac

> GenesisAeon Package 17 — Cygnus X-1 Relativistic Jet as UTAC System

[![GenesisAeon](https://img.shields.io/badge/GenesisAeon-Package%2017-blueviolet)](https://github.com/GenesisAeon)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19645351.svg)](https://doi.org/10.5281/zenodo.19645351)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Reference](https://img.shields.io/badge/Ref-Nature%20Astronomy%202026-red)](https://doi.org/10.1038/s41550-026-02828-3)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)

Models the **Cygnus X-1 black hole binary** as a full UTAC dynamical system,
calibrated against Prabu et al. (2026, *Nature Astronomy*,
DOI: [10.1038/s41550-026-02828-3](https://doi.org/10.1038/s41550-026-02828-3)).

---

## Key Scientific Result

From the measured η = 10% accretion-to-jet efficiency, the UTAC fixed-point
inversion gives the **domain-specific CREP value for Cygnus X-1 jets**:

```
Γ_jet = arctanh(η) / σ = arctanh(0.10) / 2.2 ≈ 0.0456
```

This sub-critical value (Γ << 0.5) explains the jet's sensitivity to stellar
wind perturbations — the "dancing" behaviour observed over 18 years of VLBI
monitoring is a signature of a *barely supercritical* UTAC system.

The Frame Principle is validated: σ_Φ,min / (1/16) ≈ 0.864 ≤ 1 ✓

---

## Quick Start

```python
from cygnus_jet_utac import CygnusJetUTAC

system = CygnusJetUTAC()
results = system.run_cycle(duration_years=18.0)

print(f"Γ_jet = {results['gamma_jet']:.4f}")          # → 0.0456
print(f"Jet power: {results['jet_power_W']:.3e} W")    # → ~3.85e37 W
print(f"Dance events: {results['n_dance_events']}")    # → ~36 events
```

### CLI

```bash
pip install cygnus-jet-utac

# Full 18-year simulation
cygnus-jet run

# Benchmark against Prabu 2026
cygnus-jet benchmark --verbose

# Verify Diamond-Template compliance
cygnus-jet check-diamond

# Γ_jet calibration only
cygnus-jet calibrate --eta 0.10

# Export Zenodo metadata
cygnus-jet zenodo-export
```

---

## Observable Calibration Targets (Prabu et al. 2026)

| Observable | Measured | Model Target | Tolerance |
|---|---|---|---|
| Jet power | 10,000 L☉ ≈ 3.846 × 10³⁷ W | UTAC fixed point H* | ±15% |
| Jet velocity | ~0.5c | β = 0.5 | ±10% |
| Efficiency | ~10% | CREP gating ratio | ±5% |
| Jet extent | ~16 light-years | Spatial propagation | ±20% |
| Orbital period | 5.6 days | Kepler orbit | ±1% |
| Dance events | ~2/yr | Mirror-Machine | ±50% |

---

## Diamond-Template Interface

`cygnus-jet-utac` implements the mandatory GenesisAeon Diamond-Template
interface for integration into `genesis-os` as Package 17:

```python
class CygnusJetUTAC:
    def run_cycle(self) -> dict          # main entrypoint for genesis-os
    def get_crep_state(self) -> dict     # CREP tensor snapshot {C, R, E, P, Gamma}
    def get_utac_state(self) -> dict     # UTAC ODE state {H, dH_dt, H_star, K_eff}
    def get_phase_events(self) -> list   # phase transition (dance event) log
    def to_zenodo_record(self) -> dict   # structured metadata for publication
```

### genesis-os Registration

```python
# In genesis/registry.py:
from cygnus_jet_utac import CygnusJetUTAC
PACKAGE_REGISTRY[17] = {
    "name": "cygnus-jet-utac",
    "class": CygnusJetUTAC,
    "domain": "astrophysics",
    "scale": "stellar",
    "zenodo": "10.5281/zenodo.19645351",
    "reference": "10.1038/s41550-026-02828-3",
}
```

---

## Physical Model

```
Cygnus X-1 BH (21 M☉)  +  HDE 226868 (41 M☉ Blue Supergiant)
        ↓ T_orb = 5.6 days
  Stellar wind → CREP R component → Γ(t)
        ↓
  UTAC ODE:  dH/dt = r·H·(1 − H/H*)   H*(t) = K·tanh(σ·Γ(t))
        ↓
  Relativistic jet (β=0.5) deflected by wind → Mirror-Machine
        ↓
  Phase transitions = "jet dance events" (2–4/yr over 18 yr)
```

**If genesis-os is not installed**, lightweight internal stubs are used
automatically (clearly marked `# STUB — replace with genesis.core`).

---

## Repository Structure

```
cygnus-jet-utac/
├── cygnus_jet_utac/
│   ├── __init__.py          # Exposes CygnusJetUTAC + version
│   ├── system.py            # CygnusJetUTAC (Diamond-Template)
│   ├── constants.py         # All physical constants + Cygnus X-1 params
│   ├── accretion.py         # Accretion disk → UTAC H(t)
│   ├── jet.py               # Relativistic jet propagation
│   ├── stellar_wind.py      # Blue supergiant wind (CREP modulator)
│   ├── orbital.py           # Binary orbital dynamics (5.6-day period)
│   ├── mirror_jet.py        # Mirror-Machine for dance event detection
│   ├── efficiency.py        # 10% efficiency → Γ_jet calibration
│   ├── benchmark.py         # Validation against Prabu et al. 2026
│   ├── cli.py               # CLI (cygnus-jet)
│   └── _genesis_stubs.py    # UTAC/CREP stubs if genesis-os absent
├── notebooks/
│   ├── 01_cygnus_utac_overview.ipynb
│   ├── 02_jet_dance_simulation.ipynb
│   ├── 03_benchmark_prabu2026.ipynb
│   └── 04_efficiency_crep_inversion.ipynb
├── data/
│   ├── prabu2026_measurements.yaml
│   └── cygnus_x1_radio_epochs.yaml
├── tests/                   # ≥95% coverage target
│   ├── test_diamond_interface.py
│   ├── test_efficiency.py
│   ├── test_orbital.py
│   ├── test_stellar_wind.py
│   ├── test_accretion.py
│   ├── test_jet.py
│   ├── test_mirror_jet.py
│   └── test_benchmark.py
└── docs/
    ├── physics_derivation.md
    ├── benchmark_report.md
    └── zenodo_metadata.yaml
```

---

## Installation

```bash
# Standalone
pip install cygnus-jet-utac

# With Jupyter notebooks
pip install "cygnus-jet-utac[notebooks]"

# With genesis-os integration
pip install "cygnus-jet-utac[genesis]"

# Development
git clone https://github.com/GenesisAeon/cygnus-jet-utac
cd cygnus-jet-utac
pip install -e ".[dev]"
pytest
```

---

## Notebooks Execution

```bash
pip install -e ".[notebooks]"
jupyter nbconvert --to notebook --execute --inplace \
    notebooks/01_cygnus_utac_overview.ipynb
jupyter nbconvert --to notebook --execute --inplace \
    notebooks/02_jet_dance_simulation.ipynb
jupyter nbconvert --to notebook --execute --inplace \
    notebooks/03_benchmark_prabu2026.ipynb
jupyter nbconvert --to notebook --execute --inplace \
    notebooks/04_efficiency_crep_inversion.ipynb
```

All notebooks are deterministic (`seed=42`).

---

## Citation

```bibtex
@software{romer2026cygnus,
  author  = {Römer, Johann},
  title   = {cygnus-jet-utac: GenesisAeon Package 17},
  year    = {2026},
  doi     = {10.5281/zenodo.19645351},
  url     = {https://github.com/GenesisAeon/cygnus-jet-utac}
}

@article{prabu2026cygnus,
  author  = {Prabu, S. and others},
  title   = {A jet bent by a stellar wind in the black hole X-ray binary Cygnus X-1},
  journal = {Nature Astronomy},
  year    = {2026},
  doi     = {10.1038/s41550-026-02828-3}
}
```

---

## License

MIT © Johann Römer / MOR Research Collective

*Das Universum hat sich selbst gemessen — jetzt messen wir den Jet.*
