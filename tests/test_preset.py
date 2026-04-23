"""Tests for the cygnus_jet_utac package-level interface and metadata."""

import cygnus_jet_utac
from cygnus_jet_utac import CygnusJetUTAC


def test_version_string() -> None:
    assert isinstance(cygnus_jet_utac.__version__, str)
    assert cygnus_jet_utac.__version__ == "0.1.0"


def test_zenodo_doi_present() -> None:
    assert cygnus_jet_utac.__zenodo_doi__ != ""


def test_reference_doi_present() -> None:
    assert "10.1038" in cygnus_jet_utac.__reference_doi__


def test_cygnusjetutac_instantiates() -> None:
    sys = CygnusJetUTAC()
    assert sys is not None


def test_gamma_jet_value() -> None:
    sys = CygnusJetUTAC()
    assert abs(sys.gamma_jet - 0.0456) < 0.001
