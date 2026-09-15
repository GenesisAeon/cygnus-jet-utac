# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.2] - 2026-09-15

### Fixed (documentation/metadata honesty, no numeric value change)
- **Ecosystem-wide Γ-circularity review**: `Γ_jet ≈ 0.0456` is an
  algebraic inversion of the measured 10% efficiency (Prabu et al.
  2026), not an independent calibration — it recovers the input
  exactly for ANY `sigma`. `UTAC_SIGMA_DEFAULT=2.2` (commented "ERA5
  baseline" in `constants.py`) is a shared default reused unchanged
  from the GenesisAeon climate/AMOC packages, verified byte-identical
  to `amoc-utac`'s own `UTAC_SIGMA`, not independently derived from
  Prabu et al. 2026 or any Cygnus-X-1-specific source. Further: all
  five other "Prabu 2026" benchmark targets in `benchmark.py`
  (`jet_power_W`, `jet_velocity_c`, `jet_extent_ly`,
  `orbital_period_days`, `dance_events_per_year`) were traced
  individually and are each either an echoed input constant or a free
  parameter explicitly tuned to hit that exact target — none is an
  independent prediction. The underlying UTAC ODE integration and the
  surrounding astrophysics (CAK wind law, relativistic jet kinematics,
  Keplerian orbit) are genuine, independently implemented physics —
  only the "first CREP-domain calibration" framing does not hold up.
  Removed/corrected this framing in `efficiency.py`, `benchmark.py`,
  `CITATION.cff`, `.zenodo.json`. See
  `D:\mandala\crep-utac-afet-formalism\worked_example_cygnus_jet_utac.md`
  and `FOLLOWUP_TICKETS.md` for the full analysis.

### Fixed (CI/typing only, no behavior change)
- A handful of remaining `mypy --strict` `np.ndarray`/`dict` generic
  type-parameter errors in `efficiency.py` and `benchmark.py`
  (`np.ndarray[Any, Any]`, `dict[str, Any]`), found while touching
  these files for the honesty fix above.

## [1.0.1] - 2026-09-15

### Fixed (CI only, no behavior change)
- `requires-python` was `>=3.10` while the `diamond-setup>=2.2.0`
  dependency requires `>=3.11` since 2026-06-25 — CI's `Tests (Python
  3.10)` matrix job had failed on every push since 2026-07-17 as a
  result. Bumped `requires-python` to `>=3.11`, removed the 3.10 matrix
  entry and classifier to match.
- `cygnus_jet_utac/jet.py`: `self._position`/`self._direction` had no
  explicit type annotation, so `mypy` inferred an overly narrow fixed
  shape from their `__init__` values; a later reassignment
  (`self._position = self._position + dx`) then failed `mypy --strict`'s
  shape check against a broader inferred numpy expression type. This
  had not yet surfaced in CI (no push since 2026-08-01, before the
  underlying numpy stub behavior changed) but reproduced locally with
  the current unpinned `numpy` release. Annotated both as
  `np.ndarray[Any, Any]`.
- `cygnus_jet_utac/__init__.py`'s `__version__` was never updated past
  `"0.1.0"` for the 1.0.0 release, and `tests/test_preset.py` hardcoded
  the same stale value — passing while silently testing a lie. Both
  corrected to the actual released version.

## [1.0.0] - 2026
### Added
- Initial v1.0.0 release as part of the GenesisAeon ecosystem-wide 1.0.0
  milestone.
- Standardized release tooling: `.zenodo.json`, GitHub Actions release
  workflow (`.github/workflows/release.yml`), `RELEASE_GUIDE.md`,
  `CONTRIBUTING.md`, issue/PR templates.

### Changed
- Project metadata (`pyproject.toml`) normalized: version, license,
  authors, `requires-python`, and GenesisAeon-ecosystem dependency pins
  (`genesis-os>=1.0.0`).
- License clarified as dual GPL-3.0-or-later (code) + CC BY 4.0
  (documentation and data), matching `.zenodo.json` and `README.md`;
  `LICENSE` and `pyproject.toml` previously stated MIT.
