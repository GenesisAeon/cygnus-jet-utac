# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

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
