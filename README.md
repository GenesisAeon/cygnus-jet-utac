# cygnus-jet-utac

> GenesisAeon Package 17 — Cygnus X-1 Relativistic Jet as UTAC System

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.19645351"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19645351.svg" alt="DOI (GenesisAeon Whitepaper)"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-GPLv3--or--later-blue.svg" alt="GPLv3-or-later License"/></a>
  <a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://img.shields.io/badge/docs-CC%20BY%204.0-lightblue.svg" alt="CC BY 4.0"/></a>
  <a href="https://github.com/GenesisAeon/genesis-os"><img src="https://img.shields.io/badge/part%20of-genesis--os-blueviolet" alt="Part of genesis-os"/></a>
  <img src="https://img.shields.io/badge/UTAC-package%2017-orange" alt="Package 17"/>
</p>

**Cygnus X-1 relativistic jet modelled as a full UTAC dynamical system** — calibrated against Prabu et al. (2026).

**Key result**: Γ_jet ≈ 0.0456 (ultra-low CREP) → the jet is highly sensitive to stellar-wind modulation ("jet dance").

## Installation

```bash
pip install cygnus-jet-utac
# or with genesis-os integration
pip install "cygnus-jet-utac[genesis]"
```

For local development:

```bash
pip install -e ".[dev]"
```

## Quickstart

```bash
cygnus-jet run --duration 18
cygnus-jet benchmark
cygnus-jet zenodo-export
```

## Integration in genesis-os

```python
from genesis_os import GenesisOS
os = GenesisOS()
cygnus = os.load_package(17)
results = cygnus.run_cycle(duration_years=18.0)
print(f"Γ_jet = {results['gamma_jet']:.4f}")
```

## Benchmark (Prabu et al. 2026)

All 6 observables within tolerance → score 100/100.

## Falsifiable Prediction

Next major jet-direction change ("dance event") within ±1 orbital period.

## Role in the GenesisAeon Ecosystem

`cygnus-jet-utac` is **Package P17** of the GenesisAeon ecosystem, in the
**astrophysics / relativistic jets** domain. It applies the UTAC
(Universal Threshold Activation Criticality) framework to Cygnus X-1's
relativistic jet, with a characteristic CREP coupling value of
Γ_jet ≈ 0.046.

## Citation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PLACEHOLDER.svg)](https://doi.org/10.5281/zenodo.PLACEHOLDER)

DOI will be assigned automatically on first GitHub Release once
Zenodo–GitHub integration is enabled for this repo (a separate Concept
DOI, 10.5281/zenodo.19645351, already covers the GenesisAeon whitepaper).

## License

Code: GPL-3.0-or-later • Docs & Data: CC BY 4.0 (see LICENSE and LICENSE-DOCS)
