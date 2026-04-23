"""
cygnus-jet-utac — GenesisAeon Package 17
Cygnus X-1 Relativistic Jet as a UTAC Dynamical System.

Calibrated against Prabu et al. (2026, Nature Astronomy,
DOI: 10.1038/s41550-026-02828-3).

Quick start::

    from cygnus_jet_utac import CygnusJetUTAC

    system = CygnusJetUTAC()
    results = system.run_cycle(duration_years=18.0)
    print(f"Γ_jet = {results['gamma_jet']:.4f}")
    print(f"Jet power: {results['jet_power_W']:.3e} W")
    print(f"Dance events: {len(results['phase_events'])}")
"""

from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Johann Römer"
__zenodo_doi__ = "10.5281/zenodo.19645351"
__reference_doi__ = "10.1038/s41550-026-02828-3"

from cygnus_jet_utac.system import CygnusJetUTAC

__all__ = ["CygnusJetUTAC", "__version__"]
