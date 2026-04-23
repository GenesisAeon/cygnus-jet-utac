"""CLI entrypoint for cygnus-jet-utac.

Usage::

    cygnus-jet run
    cygnus-jet run --duration 5 --dt 1800 --seed 42 --output results/
    cygnus-jet benchmark --verbose
    cygnus-jet zenodo-export --output zenodo_metadata.yaml
    cygnus-jet check-diamond
"""

from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from cygnus_jet_utac import __version__
from cygnus_jet_utac.efficiency import calibrate_gamma_jet

app = typer.Typer(
    name="cygnus-jet",
    help="GenesisAeon Package 17 — Cygnus X-1 Relativistic Jet UTAC System",
    add_completion=False,
)
console = Console()


def _make_system(dt: float = 3600.0, seed: int = 42):  # noqa: ANN201
    from cygnus_jet_utac.system import CygnusJetUTAC
    return CygnusJetUTAC(dt=dt, seed=seed)


@app.command("run")
def run_cmd(
    duration: float = typer.Option(
        18.0, "--duration", "--duration-years", help="Simulation duration (years)"
    ),
    dt: float = typer.Option(3600.0, help="Time step (s)"),
    seed: int = typer.Option(42, help="Random seed"),
    output: Path | None = typer.Option(None, help="Output directory for results"),
    plots: bool = typer.Option(False, help="Generate summary plots"),
) -> None:
    """Run the full Cygnus X-1 UTAC simulation."""
    console.print(
        Panel(
            f"[bold cyan]cygnus-jet-utac[/bold cyan] v{__version__}\n"
            f"GenesisAeon Package 17 · Cygnus X-1 UTAC Simulation\n"
            f"Duration: {duration} yr · dt: {dt:.0f} s · seed: {seed}",
            title="[bold]Starting Simulation[/bold]",
            border_style="blue",
        )
    )

    system = _make_system(dt=dt, seed=seed)

    with console.status(
        f"[bold green]Simulating {duration} years ({int(duration*365.25*86400/dt):,} steps)...",
        spinner="dots",
    ):
        results = system.run_cycle(duration_years=duration)

    # Summary table
    tbl = Table(title="Simulation Results", show_header=True, header_style="bold magenta")
    tbl.add_column("Observable", style="cyan")
    tbl.add_column("Value", justify="right")
    tbl.add_column("Unit")

    tbl.add_row("Γ_jet (calibrated)", f"{results['gamma_jet']:.4f}", "dimensionless")
    tbl.add_row("σ_Φ,min / (1/16)", f"{results['sigma_phi_ratio']:.4f}", "dimensionless")
    tbl.add_row(
        "Frame Principle",
        "✓ satisfied" if results["frame_principle_satisfied"] else "✗ violated",
        "",
    )
    tbl.add_row("Jet power", f"{results['jet_power_W']:.3e}", "W")
    tbl.add_row("Jet power", f"{results['jet_power_Lsun']:.0f}", "L☉")
    tbl.add_row("Jet velocity β", f"{results['jet_velocity_c']:.2f}", "c")
    tbl.add_row("Jet extent", f"{results['jet_extent_ly']:.2f}", "ly")
    tbl.add_row("Orbital period", f"{results['orbital_period_days']:.2f}", "days")
    tbl.add_row("Dance events", str(results["n_dance_events"]), "total")
    tbl.add_row("Dance rate", f"{results['dance_events_per_year']:.2f}", "events/yr")
    tbl.add_row("Genesis-OS", "available" if results["genesis_available"] else "stub mode", "")

    console.print(tbl)

    if output:
        output.mkdir(parents=True, exist_ok=True)
        out_file = output / "results.json"
        # Serialise (drop large list histories for the summary file)
        slim = {k: v for k, v in results.items()
                if not isinstance(v, list) or k == "phase_events"}
        out_file.write_text(json.dumps(slim, indent=2, default=str))
        console.print(f"[green]Results saved to {out_file}")

        if plots:
            system.plot_summary(save_path=str(output / "summary.png"))
            system.plot_jet_trajectory(save_path=str(output / "jet_trajectory.png"))
            system.plot_crep_evolution(save_path=str(output / "crep_evolution.png"))
            system.plot_phase_events(save_path=str(output / "dance_events.png"))
            console.print(f"[green]Plots saved to {output}/")


@app.command("benchmark")
def benchmark_cmd(
    verbose: bool = typer.Option(False, help="Print full benchmark report"),
    dt: float = typer.Option(3600.0, help="Time step (s)"),
    seed: int = typer.Option(42, help="Random seed"),
    duration_years: float = typer.Option(18.0, help="Simulation duration (years)"),
) -> None:
    """Validate simulation against Prabu et al. 2026 measurements."""
    from cygnus_jet_utac.benchmark import run_benchmark

    console.print("[bold cyan]Running benchmark against Prabu et al. 2026...[/bold cyan]")
    system = _make_system(dt=dt, seed=seed)

    with console.status(f"[bold green]Simulating {duration_years} years...", spinner="dots"):
        system.run_cycle(duration_years=duration_years)

    bm = run_benchmark(system)

    if verbose:
        console.print(bm["report_str"])
    else:
        score_color = "green" if bm["score"] >= 80 else ("yellow" if bm["score"] >= 60 else "red")
        console.print(
            Panel(
                f"[bold {score_color}]Score: {bm['score']:.1f}/100[/bold {score_color}]\n"
                f"Passed: {bm['n_passed']}/{bm['n_total']} observables\n"
                f"Γ_jet = {bm['gamma_jet']:.4f}\n"
                f"Frame Principle: {'✓' if bm['sigma_phi_satisfied'] else '✗'}",
                title="Benchmark Result",
                border_style=score_color,
            )
        )

    if bm["score"] < 80:
        console.print("[yellow]Score below 80/100 — consider tuning parameters.[/yellow]")


@app.command("zenodo-export")
def zenodo_export_cmd(
    output: Path = typer.Option(
        Path("zenodo_metadata.json"), help="Output JSON file path"
    ),
    dt: float = typer.Option(3600.0, help="Time step (s)"),
    seed: int = typer.Option(42, help="Random seed"),
    duration_years: float = typer.Option(18.0, help="Simulation duration (years)"),
) -> None:
    """Export Zenodo-ready metadata for this simulation run."""
    console.print("[bold cyan]Generating Zenodo metadata...[/bold cyan]")
    system = _make_system(dt=dt, seed=seed)

    with console.status(f"[bold green]Running simulation ({duration_years} yr)...", spinner="dots"):
        system.run_cycle(duration_years=duration_years)

    record = system.to_zenodo_record()
    output.write_text(json.dumps(record, indent=2, default=str))
    console.print(f"[green]Zenodo record saved to {output}")
    console.print(f"  DOI: {record.get('doi')}")


@app.command("check-diamond")
def check_diamond_cmd() -> None:
    """Verify Diamond-Template interface compliance."""
    from cygnus_jet_utac.system import CygnusJetUTAC

    console.print("[bold cyan]Checking Diamond-Template compliance...[/bold cyan]")
    system = CygnusJetUTAC(dt=86400.0)

    required_methods = [
        "run_cycle", "get_crep_state", "get_utac_state",
        "get_phase_events", "to_zenodo_record",
    ]
    all_ok = True
    for method in required_methods:
        ok = callable(getattr(system, method, None))
        color = "green" if ok else "red"
        console.print(f"  [{color}]{'✓' if ok else '✗'}[/{color}]  {method}()")
        all_ok = all_ok and ok

    # Quick functional checks
    crep = system.get_crep_state()
    utac = system.get_utac_state()
    events = system.get_phase_events()
    record = system.to_zenodo_record()

    checks = [
        ("CREP keys present", all(k in crep for k in ["C", "R", "E", "P", "Gamma"])),
        ("UTAC keys present", all(k in utac for k in ["H", "dH_dt", "H_star", "K_eff"])),
        ("phase_events is list", isinstance(events, list)),
        ("Zenodo DOI present", "doi" in record),
        ("Γ_jet ≈ 0.0456", abs(system.gamma_jet - 0.0456) < 0.001),
    ]
    for desc, ok in checks:
        color = "green" if ok else "red"
        console.print(f"  [{color}]{'✓' if ok else '✗'}[/{color}]  {desc}")
        all_ok = all_ok and ok

    status = "SATISFIED ✓" if all_ok else "FAILED ✗"
    color = "green" if all_ok else "red"
    console.print(f"\n[bold {color}]Diamond-Template contract: {status}[/bold {color}]")


@app.command("calibrate")
def calibrate_cmd(
    eta: float = typer.Option(0.10, help="Accretion-to-jet efficiency"),
    sigma: float = typer.Option(2.2, help="UTAC CREP coupling σ"),
) -> None:
    """Print Γ_jet calibration from efficiency inversion."""
    result = calibrate_gamma_jet(eta=eta, sigma=sigma, verbose=True)
    console.print(result["report"])


def main() -> None:
    """Main CLI entrypoint."""
    app()


if __name__ == "__main__":
    main()
