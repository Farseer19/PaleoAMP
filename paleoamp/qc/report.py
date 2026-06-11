"""Generate QC reports from DamageProfile results."""

from __future__ import annotations

import json
from pathlib import Path

from rich.console import Console
from rich.table import Table
from rich import box

from paleoamp.qc.adna import DamageProfile


_console = Console()


def print_summary_table(profiles: list[DamageProfile]) -> None:
    """Print a Rich table summarising QC results for all assessed files."""
    table = Table(
        title="PaleoAMP — Ancient DNA Quality Report",
        box=box.ROUNDED,
        highlight=True,
    )
    table.add_column("File", style="cyan", no_wrap=False, max_width=40)
    table.add_column("Reads", justify="right")
    table.add_column("Passed\nLen filter", justify="right")
    table.add_column("Mean\nPhred", justify="right")
    table.add_column("Median\nLength (bp)", justify="right")
    table.add_column("5′ C→T\n(pos 1)", justify="right")
    table.add_column("3′ G→A\n(pos 1)", justify="right")
    table.add_column("Verdict", justify="center")

    for p in profiles:
        ct = p.ct_rates_5prime[0] if p.ct_rates_5prime else 0.0
        ga = p.ga_rates_3prime[0] if p.ga_rates_3prime else 0.0
        verdict = "[green]PASS[/]" if p.pass_verdict else "[red]FAIL[/]"
        frac = (
            f"{p.passed_length_filter / p.total_reads:.1%}"
            if p.total_reads > 0 else "—"
        )
        table.add_row(
            Path(p.path).name,
            f"{p.total_reads:,}",
            frac,
            f"{p.mean_phred:.1f}" if p.mean_phred else "—",
            f"{p.median_read_length:.0f}" if p.median_read_length else "—",
            f"{ct:.3f}",
            f"{ga:.3f}",
            verdict,
        )

    _console.print(table)

    failed = [p for p in profiles if not p.pass_verdict]
    if failed:
        _console.print("\n[bold red]Failure reasons:[/]")
        for p in failed:
            _console.print(f"  [cyan]{Path(p.path).name}[/]")
            for reason in p.fail_reasons:
                _console.print(f"    • {reason}")


def print_damage_plot(profile: DamageProfile, positions: int = 10) -> None:
    """Print ASCII misincorporation plot for a single profile."""
    _console.rule(f"[bold]{Path(profile.path).name}[/] — Damage Profile")

    ct = profile.ct_rates_5prime[:positions]
    ga = profile.ga_rates_3prime[:positions]

    _console.print("\n5′ end  C→T misincorporation rate (pos 1 = leftmost base):")
    _print_bar_chart(ct, label_prefix="5′ pos")

    _console.print("\n3′ end  G→A misincorporation rate (pos 1 = rightmost base):")
    _print_bar_chart(ga, label_prefix="3′ pos")
    _console.print()


def _print_bar_chart(rates: list[float], label_prefix: str) -> None:
    max_width = 40
    for i, rate in enumerate(rates, start=1):
        bar_len = int(rate * max_width / max(rates + [0.001]))
        bar = "█" * bar_len
        color = "green" if rate >= 0.05 else "yellow" if rate >= 0.02 else "red"
        _console.print(
            f"  {label_prefix} {i:>2}: [{color}]{bar:<{max_width}}[/] {rate:.4f}"
        )


def save_json_report(profiles: list[DamageProfile], output_path: Path) -> None:
    """Save full QC results to a JSON file for downstream use."""
    records = []
    for p in profiles:
        records.append({
            "path": p.path,
            "total_reads": p.total_reads,
            "passed_length_filter": p.passed_length_filter,
            "mean_phred": round(p.mean_phred, 2),
            "mean_read_length": round(p.mean_read_length, 1),
            "median_read_length": round(p.median_read_length, 1),
            "ct_rates_5prime": [round(r, 5) for r in p.ct_rates_5prime],
            "ga_rates_3prime": [round(r, 5) for r in p.ga_rates_3prime],
            "pass_verdict": p.pass_verdict,
            "fail_reasons": p.fail_reasons,
        })

    with open(output_path, "w") as fh:
        json.dump(records, fh, indent=2)
    _console.print(f"\n[green]JSON report saved → {output_path}[/]")


def save_tsv_report(profiles: list[DamageProfile], output_path: Path) -> None:
    """Save a flat TSV summary suitable for downstream filtering."""
    header = (
        "file\ttotal_reads\tpassed_length_filter\tmean_phred\t"
        "mean_read_length\tmedian_read_length\t"
        "ct_rate_5p_pos1\tga_rate_3p_pos1\tpass_verdict\tfail_reasons\n"
    )
    with open(output_path, "w") as fh:
        fh.write(header)
        for p in profiles:
            ct = p.ct_rates_5prime[0] if p.ct_rates_5prime else 0.0
            ga = p.ga_rates_3prime[0] if p.ga_rates_3prime else 0.0
            reasons = "; ".join(p.fail_reasons) or ""
            fh.write(
                f"{p.path}\t{p.total_reads}\t{p.passed_length_filter}\t"
                f"{p.mean_phred:.2f}\t{p.mean_read_length:.1f}\t"
                f"{p.median_read_length:.1f}\t"
                f"{ct:.5f}\t{ga:.5f}\t{p.pass_verdict}\t{reasons}\n"
            )
    _console.print(f"[green]TSV report saved → {output_path}[/]")
