"""Tests for the cygnus-jet-utac CLI commands."""

import json

from typer.testing import CliRunner

from cygnus_jet_utac.cli import app

runner = CliRunner()


def test_calibrate_command() -> None:
    result = runner.invoke(app, ["calibrate"])
    assert result.exit_code == 0
    assert "0.04" in result.output or "gamma" in result.output.lower()


def test_check_diamond_command() -> None:
    result = runner.invoke(app, ["check-diamond"])
    assert result.exit_code == 0
    assert "diamond" in result.output.lower() or "pass" in result.output.lower()


def test_run_command_short() -> None:
    result = runner.invoke(app, ["run", "--duration-years", "0.1", "--dt", "86400"])
    assert result.exit_code == 0, result.output


def test_benchmark_command_short() -> None:
    result = runner.invoke(
        app, ["benchmark", "--duration-years", "0.1", "--dt", "86400"]
    )
    assert result.exit_code == 0, result.output
    assert "score" in result.output.lower() or "benchmark" in result.output.lower()


def test_zenodo_export_command(tmp_path) -> None:
    out_file = tmp_path / "zenodo.json"
    result = runner.invoke(
        app,
        [
            "zenodo-export",
            "--duration-years", "0.1",
            "--dt", "86400",
            "--output", str(out_file),
        ],
    )
    assert result.exit_code == 0, result.output
    assert out_file.exists()
    data = json.loads(out_file.read_text())
    assert "title" in data or "version" in data or "doi" in data
