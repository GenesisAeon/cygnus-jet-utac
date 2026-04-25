# cygnus-jet-utac

> GenesisAeon Package 17 — Cygnus X-1 Relativistic Jet as UTAC System

[![GenesisAeon](https://img.shields.io/badge/GenesisAeon-Package%2017-blueviolet)](https://github.com/GenesisAeon)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19645351.svg)](https://doi.org/10.5281/zenodo.19645351)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Reference](https://img.shields.io/badge/Ref-Nature%20Astronomy%202026-red)](https://doi.org/10.1038/s41550-026-02828-3)

**Cygnus X-1 relativistic jet modelled as a full UTAC dynamical system** — calibrated against Prabu et al. (2026).

**Key result**: Γ_jet ≈ 0.0456 (ultra-low CREP) → the jet is highly sensitive to stellar-wind modulation ("jet dance").

## Installation

```bash
pip install -e ".[dev]"
# or with genesis-os integration
pip install -e ".[genesis]"
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

## License

Code: MIT • Docs & Data: CC BY 4.0
