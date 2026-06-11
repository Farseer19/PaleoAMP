"""Command-line interface for PaleoAMP."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import yaml
from rich.console import Console

_console = Console()


def _load_config(config_path: Path | None) -> dict:
    if config_path is None:
        default = Path(__file__).parent.parent / "config" / "defaults.yaml"
        config_path = default if default.exists() else None

    if config_path and Path(config_path).exists():
        with open(config_path) as fh:
            return yaml.safe_load(fh) or {}
    return {}


@click.group()
@click.version_option()
def main():
    """PaleoAMP — discover novel AMP genes in ancient/extinct microbial samples."""


# ---------------------------------------------------------------------------
# paleoamp db
# ---------------------------------------------------------------------------

@main.group("db")
def db_group():
    """Manage AMP reference databases."""


@db_group.command("download")
@click.option(
    "--output-dir", "-o",
    default="./data/amp_databases",
    show_default=True,
    help="Directory to save downloaded database files.",
)
@click.option(
    "--databases", "-d",
    multiple=True,
    default=["dramp", "apd3", "dbamp"],
    show_default=True,
    help="Which databases to download. Repeat flag for multiple.",
)
@click.option("--force", is_flag=True, help="Re-download even if file exists.")
def db_download(output_dir: str, databases: tuple, force: bool):
    """Download AMP reference databases (DRAMP, APD3, dbAMP)."""
    from paleoamp.data.amp_db import download_amp_databases

    paths = download_amp_databases(
        output_dir=Path(output_dir),
        databases=list(databases),
        force=force,
    )
    _console.print(f"\n[bold green]{len(paths)} database(s) ready.[/]")


@db_group.command("merge")
@click.option(
    "--input-dir", "-i",
    default="./data/amp_databases",
    show_default=True,
    help="Directory containing downloaded database FASTA files.",
)
@click.option(
    "--output", "-o",
    default="./data/amp_databases/merged_amp_db.fasta",
    show_default=True,
)
@click.option("--no-dedup", is_flag=True, help="Skip deduplication.")
def db_merge(input_dir: str, output: str, no_dedup: bool):
    """Merge all downloaded AMP databases into a single non-redundant FASTA."""
    from paleoamp.data.amp_db import merge_amp_databases

    in_dir = Path(input_dir)
    fasta_files = {
        p.stem: p
        for p in in_dir.glob("*.fasta")
        if "merged" not in p.name
    }
    if not fasta_files:
        _console.print(f"[red]No FASTA files found in {in_dir}[/]")
        sys.exit(1)

    merge_amp_databases(
        db_paths=fasta_files,
        output_path=Path(output),
        deduplicate=not no_dedup,
    )


# ---------------------------------------------------------------------------
# paleoamp sra
# ---------------------------------------------------------------------------

@main.group("sra")
def sra_group():
    """Search and download ancient metagenome data from NCBI SRA."""


@sra_group.command("search")
@click.argument("query")
@click.option("--max-results", "-n", default=50, show_default=True)
@click.option(
    "--output", "-o",
    default=None,
    help="Save results to this TSV file.",
)
@click.option(
    "--ancient-filter/--no-ancient-filter",
    default=True,
    show_default=True,
    help="Apply heuristic keyword filter for ancient DNA indicators.",
)
@click.option(
    "--microbial-only",
    is_flag=True,
    default=False,
    help=(
        "Keep only results with METAGENOMIC library source and microbial "
        "context keywords. Excludes plant/animal eDNA studies."
    ),
)
@click.option(
    "--min-reads",
    default=None,
    type=int,
    metavar="N",
    help=(
        "Discard runs with fewer than N reads (total_spots). "
        "Recommended ≥500000 for a workable aDNA assembly."
    ),
)
@click.option(
    "--min-bases",
    default=None,
    type=int,
    metavar="BP",
    help="Discard runs with fewer than BP total bases.",
)
@click.option(
    "--min-size-mb",
    default=None,
    type=float,
    metavar="MB",
    help="Discard runs whose reported size is below MB megabytes.",
)
def sra_search(
    query: str,
    max_results: int,
    output: str | None,
    ancient_filter: bool,
    microbial_only: bool,
    min_reads: int | None,
    min_bases: int | None,
    min_size_mb: float | None,
):
    """Search SRA for ancient metagenome datasets matching QUERY."""
    from paleoamp.data.sra import (
        search_ancient_metagenomes,
        filter_for_ancient_indicators,
        filter_for_microbial,
        filter_for_assembly_size,
    )

    _console.print(f'[bold]Searching SRA:[/] "{query}"')
    df = search_ancient_metagenomes(query, max_results=max_results)

    if df.empty:
        _console.print("[yellow]No results found.[/]")
        return

    if ancient_filter:
        before = len(df)
        df = filter_for_ancient_indicators(df)
        _console.print(f"Ancient filter:   {len(df):>4} / {before} results retained.")

    if microbial_only:
        before = len(df)
        df = filter_for_microbial(df)
        _console.print(f"Microbial filter: {len(df):>4} / {before} results retained.")

    if min_reads is not None or min_bases is not None or min_size_mb is not None:
        before = len(df)
        df = filter_for_assembly_size(df, min_spots=min_reads, min_bases=min_bases, min_size_mb=min_size_mb)
        _console.print(f"Size filter:      {len(df):>4} / {before} results retained.")

    if df.empty:
        _console.print("[yellow]No results after filtering.[/]")
        return

    _console.print(df.to_string(index=False))

    if output:
        df.to_csv(output, sep="\t", index=False)
        _console.print(f"\n[green]Results saved → {output}[/]")


@sra_group.command("merge-candidates")
@click.argument("inputs", nargs=-1, required=True, metavar="TSV [TSV ...]")
@click.option(
    "--output", "-o",
    default="./candidates_merged.tsv",
    show_default=True,
    help="Output path for the merged TSV.",
)
@click.option(
    "--sort-by",
    default="total_spots",
    show_default=True,
    help="Column to sort results by (descending).",
)
def sra_merge_candidates(inputs: tuple, output: str, sort_by: str):
    """
    Merge and deduplicate candidate TSV files from multiple searches.

    Accepts any number of TSV files (or glob patterns) produced by
    'paleoamp sra search'. Deduplicates on run_accession and adds a
    'found_in' column listing every source file each run appeared in.

    Example:

        paleoamp sra merge-candidates results/*.tsv -o merged.tsv
    """
    import glob

    import pandas as pd

    # Expand any glob patterns and collect all paths
    all_paths: list[Path] = []
    for pattern in inputs:
        expanded = glob.glob(pattern)
        if expanded:
            all_paths.extend(Path(p) for p in sorted(expanded))
        else:
            all_paths.append(Path(pattern))

    missing = [p for p in all_paths if not p.exists()]
    if missing:
        for p in missing:
            _console.print(f"[red]Not found: {p}[/]")
        sys.exit(1)

    if not all_paths:
        _console.print("[red]No input files found.[/]")
        sys.exit(1)

    # Load each file, tag with its source filename
    frames: list[pd.DataFrame] = []
    total_rows = 0
    for path in all_paths:
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception as exc:
            _console.print(f"[yellow]Skipping {path.name}: {exc}[/]")
            continue
        if df.empty or "run_accession" not in df.columns:
            _console.print(f"[yellow]Skipping {path.name}: no run_accession column[/]")
            continue
        df["_source"] = path.name
        total_rows += len(df)
        frames.append(df)
        _console.print(f"  loaded [cyan]{path.name}[/]: {len(df):,} rows")

    if not frames:
        _console.print("[red]No valid input files.[/]")
        sys.exit(1)

    combined = pd.concat(frames, ignore_index=True)

    # Build found_in: comma-separated list of source files per run_accession
    found_in = (
        combined.groupby("run_accession")["_source"]
        .apply(lambda s: ", ".join(sorted(s.unique())))
        .rename("found_in")
    )

    # Deduplicate: keep first occurrence of each run_accession
    deduped = (
        combined.drop(columns=["_source"])
        .drop_duplicates(subset=["run_accession"], keep="first")
        .merge(found_in, on="run_accession", how="left")
    )

    # Sort
    if sort_by in deduped.columns:
        deduped = deduped.sort_values(sort_by, ascending=False, na_position="last")
    else:
        _console.print(f"[yellow]Sort column '{sort_by}' not found; skipping sort.[/]")

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    deduped.to_csv(out_path, sep="\t", index=False)

    duplicates_removed = total_rows - len(deduped)
    in_multiple = (found_in.str.contains(",")).sum()

    _console.print(
        f"\n[bold green]Merged {len(all_paths)} file(s):[/] "
        f"{total_rows:,} total rows → "
        f"[bold]{len(deduped):,} unique runs[/] "
        f"({duplicates_removed:,} duplicates removed)"
    )
    if in_multiple:
        _console.print(
            f"  {in_multiple:,} run(s) appeared in more than one search."
        )
    _console.print(f"  Saved → [cyan]{out_path}[/]")


@sra_group.command("download")
@click.argument("accessions", nargs=-1, required=True)
@click.option(
    "--output-dir", "-o",
    default="./data/reads",
    show_default=True,
)
@click.option("--threads", "-t", default=4, show_default=True)
@click.option(
    "--via-ena/--via-sratools",
    default=True,
    show_default=True,
    help="Download via ENA HTTP (no SRA Toolkit needed) or fasterq-dump.",
)
def sra_download(accessions: tuple, output_dir: str, threads: int, via_ena: bool):
    """Download FASTQ files for one or more SRA accessions."""
    if via_ena:
        from paleoamp.data.sra import download_via_ena, get_run_metadata
        _console.print("Fetching ENA URLs from metadata ...")
        meta = get_run_metadata(list(accessions))
        results = download_via_ena(
            accessions=list(accessions),
            output_dir=Path(output_dir),
            metadata_df=meta if not meta.empty else None,
        )
    else:
        from paleoamp.data.sra import download_runs
        results = download_runs(
            accessions=list(accessions),
            output_dir=Path(output_dir),
            threads=threads,
        )
    _console.print(f"\n[bold green]{len(results)}/{len(accessions)} run(s) downloaded.[/]")
    for acc, path in results.items():
        _console.print(f"  {acc} → {path}")


# ---------------------------------------------------------------------------
# paleoamp assemble
# ---------------------------------------------------------------------------

@main.group("assemble")
def assemble_group():
    """Assemble QC-passed ancient metagenome reads into contigs."""


@assemble_group.command("run")
@click.option(
    "--reads-dir", "-r",
    default="./data/reads",
    show_default=True,
    help="Directory containing per-sample FASTQ subdirectories.",
)
@click.option(
    "--qc-report",
    default="./results/qc/qc_report.json",
    show_default=True,
    type=click.Path(exists=True),
    help="QC JSON report produced by 'paleoamp qc assess'.",
)
@click.option(
    "--output-dir", "-o",
    default="./results/assembly",
    show_default=True,
    help="Directory to write per-sample MEGAHIT output.",
)
@click.option("--threads", "-t", default=4, show_default=True)
@click.option(
    "--min-contig-len",
    default=200,
    show_default=True,
    help="Discard contigs shorter than this (bp).",
)
def assemble_run(
    reads_dir: str,
    qc_report: str,
    output_dir: str,
    threads: int,
    min_contig_len: int,
):
    """
    Run MEGAHIT on each sample that passed aDNA QC.

    Reads the QC report to determine which samples to assemble and
    automatically skips any that failed quality assessment.
    """
    from paleoamp.assembly.assemble import assemble_all
    from rich.table import Table
    from rich import box

    _console.print(f"Reading QC verdicts from [cyan]{qc_report}[/]")

    results = assemble_all(
        reads_dir=Path(reads_dir),
        output_dir=Path(output_dir),
        qc_report=Path(qc_report),
        threads=threads,
        min_contig_len=min_contig_len,
    )

    table = Table(
        title="PaleoAMP — Assembly Results",
        box=box.ROUNDED,
        highlight=True,
    )
    table.add_column("Sample", style="cyan")
    table.add_column("Contigs", justify="right")
    table.add_column("Total (bp)", justify="right")
    table.add_column("N50 (bp)", justify="right")
    table.add_column("Status", justify="center")

    assembled = 0
    for r in results:
        if r.skipped:
            status = f"[yellow]SKIPPED — {r.skip_reason}[/]"
            table.add_row(r.sample_id, "—", "—", "—", status)
        else:
            assembled += 1
            status = "[green]OK[/]"
            table.add_row(
                r.sample_id,
                f"{r.n_contigs:,}",
                f"{r.total_length_bp:,}",
                f"{r.n50:,}",
                status,
            )

    _console.print(table)
    _console.print(
        f"\n[bold]{assembled}/{len(results)} sample(s) assembled.[/] "
        f"Contigs in [cyan]{output_dir}[/]"
    )
    if assembled == 0:
        sys.exit(1)


# ---------------------------------------------------------------------------
# paleoamp predict
# ---------------------------------------------------------------------------

@main.group("predict")
def predict_group():
    """Predict ORFs from assembled contigs."""


@predict_group.command("orfs")
@click.option(
    "--assembly-dir", "-a",
    default="./results/assembly",
    show_default=True,
    help="Directory of per-sample MEGAHIT output (from 'paleoamp assemble run').",
)
@click.option(
    "--output-dir", "-o",
    default="./results/orfs",
    show_default=True,
    help="Directory to write per-sample .faa and .gff files.",
)
def predict_orfs_cmd(assembly_dir: str, output_dir: str):
    """
    Run Prodigal (metagenome mode) on each assembled sample.

    Writes a protein FASTA (.faa) and GFF coordinate file (.gff) per sample.
    """
    from paleoamp.predict.orf import predict_all
    from rich.table import Table
    from rich import box

    results = predict_all(
        assembly_dir=Path(assembly_dir),
        output_dir=Path(output_dir),
    )

    table = Table(
        title="PaleoAMP — ORF Prediction Results",
        box=box.ROUNDED,
        highlight=True,
    )
    table.add_column("Sample", style="cyan")
    table.add_column("ORFs predicted", justify="right")
    table.add_column("Protein FASTA", style="dim")
    table.add_column("Status", justify="center")

    predicted = 0
    for r in results:
        if r.skipped:
            table.add_row(r.sample_id, "—", "—", f"[yellow]SKIPPED — {r.skip_reason}[/]")
        else:
            predicted += 1
            table.add_row(
                r.sample_id,
                f"{r.n_orfs:,}",
                str(r.protein_fasta),
                "[green]OK[/]",
            )

    _console.print(table)
    _console.print(
        f"\n[bold]{predicted}/{len(results)} sample(s) processed.[/] "
        f"ORFs in [cyan]{output_dir}[/]"
    )
    if predicted == 0:
        sys.exit(1)


@predict_group.command("hydrophobic")
@click.option(
    "--orfs-dir", "-i",
    default="./results/orfs",
    show_default=True,
    help="Directory of per-sample .faa files (from 'paleoamp predict orfs').",
)
@click.option(
    "--output-dir", "-o",
    default="./results/orfs_hydrophobic",
    show_default=True,
    help="Directory to write filtered .faa files.",
)
@click.option(
    "--min-kd",
    default=0.0,
    show_default=True,
    type=float,
    help="Minimum mean Kyte-Doolittle hydrophobicity score.",
)
@click.option(
    "--min-hydrophobic-fraction",
    default=0.30,
    show_default=True,
    type=float,
    help="Minimum fraction of hydrophobic residues (ACFILMVWY).",
)
def predict_hydrophobic_cmd(
    orfs_dir: str,
    output_dir: str,
    min_kd: float,
    min_hydrophobic_fraction: float,
):
    """
    Filter predicted ORFs by hydrophobicity as a pre-screen for AMP candidates.

    Retains ORFs whose mean Kyte-Doolittle score and hydrophobic residue
    fraction both meet the given thresholds. Writes one filtered .faa per
    sample to OUTPUT_DIR for downstream BLAST classification.
    """
    from paleoamp.predict.hydrophobic import filter_by_hydrophobicity
    from rich.table import Table
    from rich import box

    orfs_path = Path(orfs_dir)
    out_path = Path(output_dir)

    sample_dirs = sorted(d for d in orfs_path.iterdir() if d.is_dir())
    if not sample_dirs:
        _console.print(f"[red]No sample directories found in {orfs_path}[/]")
        sys.exit(1)

    table = Table(
        title="PaleoAMP — Hydrophobicity Filter",
        box=box.ROUNDED,
        highlight=True,
    )
    table.add_column("Sample", style="cyan")
    table.add_column("Input ORFs", justify="right")
    table.add_column("Passed", justify="right")
    table.add_column("Retained %", justify="right")

    total_in = total_out = 0
    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        faa = sample_dir / f"{sample_id}.faa"
        if not faa.exists():
            continue
        out_faa = out_path / sample_id / f"{sample_id}.faa"
        n_in, n_out = filter_by_hydrophobicity(
            faa, out_faa,
            min_kd=min_kd,
            min_hydrophobic_fraction=min_hydrophobic_fraction,
        )
        total_in += n_in
        total_out += n_out
        pct = f"{n_out / n_in * 100:.1f}%" if n_in else "—"
        table.add_row(sample_id, f"{n_in:,}", f"{n_out:,}", pct)

    _console.print(table)
    overall = f"{total_out / total_in * 100:.1f}%" if total_in else "—"
    _console.print(
        f"\n[bold]{total_out:,} / {total_in:,} ORFs retained ({overall})[/] "
        f"→ [cyan]{output_dir}[/]"
    )


# ---------------------------------------------------------------------------
# paleoamp ml
# ---------------------------------------------------------------------------

def _latest_sample_orfs(orfs_dir: Path) -> Path:
    """
    Return the .faa path for the sample whose FAA file was most recently
    created under *orfs_dir*.  Uses the FAA mtime (set when Prodigal writes
    it) rather than the directory mtime, which changes on any access.
    """
    faa_files = [
        d / f"{d.name}.faa"
        for d in orfs_dir.iterdir()
        if d.is_dir() and (d / f"{d.name}.faa").exists()
    ]
    if not faa_files:
        raise FileNotFoundError(f"No .faa files found under {orfs_dir}")
    return max(faa_files, key=lambda f: f.stat().st_mtime)


@main.group("ml")
def ml_group():
    """Train and run the ESM-2 + MLP AMP classifier."""


@ml_group.command("collect")
@click.option(
    "--output-dir", "-o",
    default="./data/ml",
    show_default=True,
    help="Directory to download raw source files and write the merged dataset.",
)
@click.option(
    "--sources", "-s",
    multiple=True,
    default=[],
    help="Specific sources to download (default: all). Repeat for multiple.",
)
@click.option("--force", is_flag=True, help="Re-download even if files exist.")
@click.option("--min-len", default=5, show_default=True)
@click.option("--max-len", default=200, show_default=True)
def ml_collect(output_dir: str, sources: tuple, force: bool, min_len: int, max_len: int):
    """
    Download positive and negative AMP training sequences and build the
    merged, filtered, deduplicated dataset CSV.

    Sources — positives: apd6, dramp, dbaasp
               negatives: amplify_neg, uniprot_non_amp
    """
    from paleoamp.ml.collect import download_sources, build_dataset, save_dataset

    raw_dir = Path(output_dir) / "raw"
    _console.print("[bold]Downloading AMP training sources ...[/]")
    download_sources(raw_dir, list(sources) or None, force=force)

    _console.print("\n[bold]Building merged dataset ...[/]")
    df = build_dataset(raw_dir, min_len=min_len, max_len=max_len)

    out_csv = Path(output_dir) / "dataset.csv"
    save_dataset(df, out_csv)
    _console.print(f"\n[bold green]Dataset ready:[/] {len(df):,} sequences → [cyan]{out_csv}[/]")


@ml_group.command("cluster")
@click.option(
    "--dataset", "-i",
    default="./data/ml/dataset.csv",
    show_default=True,
    type=click.Path(exists=True),
    help="Dataset CSV from 'paleoamp ml collect'.",
)
@click.option(
    "--output-dir", "-o",
    default="./data/ml",
    show_default=True,
)
@click.option("--min-seq-id", default=0.40, show_default=True, type=float,
              help="Maximum identity threshold for clustering (0.0–1.0).")
@click.option("--test-fraction", default=0.20, show_default=True, type=float)
@click.option("--threads", "-t", default=4, show_default=True)
@click.option("--seed", default=42, show_default=True)
def ml_cluster(
    dataset: str,
    output_dir: str,
    min_seq_id: float,
    test_fraction: float,
    threads: int,
    seed: int,
):
    """
    Cluster the dataset at --min-seq-id identity then assign whole clusters
    to train/test splits (cluster-based split prevents data leakage).
    """
    from paleoamp.ml.collect import load_dataset, save_dataset
    from paleoamp.ml.cluster import cluster_sequences, cluster_split

    _console.print("[bold]Loading dataset ...[/]")
    df = load_dataset(Path(dataset))

    work_dir = Path(output_dir) / "cluster_work"
    _console.print(f"[bold]Clustering {len(df):,} sequences at {min_seq_id:.0%} identity ...[/]")
    cluster_ids = cluster_sequences(df, work_dir, min_seq_id=min_seq_id, threads=threads)

    df_split = cluster_split(df, cluster_ids, test_fraction=test_fraction, random_seed=seed)
    out = Path(output_dir) / "dataset_split.csv"
    save_dataset(df_split, out)
    _console.print(f"[bold green]Split dataset saved → [cyan]{out}[/]")


@ml_group.command("embed")
@click.option(
    "--dataset", "-i",
    default="./data/ml/dataset_split.csv",
    show_default=True,
    type=click.Path(exists=True),
    help="Split dataset CSV from 'paleoamp ml cluster'.",
)
@click.option(
    "--output-dir", "-o",
    default="./data/ml/embeddings",
    show_default=True,
)
@click.option("--batch-size", default=64, show_default=True)
@click.option("--device", default="auto", show_default=True,
              help="'auto', 'cpu', or 'cuda'.")
def ml_embed(dataset: str, output_dir: str, batch_size: int, device: str):
    """
    Generate frozen ESM-2 embeddings for train and test splits.

    Embeddings are saved as .npy files and cached — re-running is instant.
    """
    from paleoamp.ml.collect import load_dataset
    from paleoamp.ml.embed import embed_dataframe

    df = load_dataset(Path(dataset))
    out = Path(output_dir)
    for split in ("train", "test"):
        n = (df["split"] == split).sum()
        _console.print(f"\n[bold]Embedding {split} split ({n:,} sequences) ...[/]")
        embed_dataframe(df, out, batch_size=batch_size, device=device, split=split)
    _console.print(f"\n[bold green]Embeddings saved → [cyan]{out}[/]")


@ml_group.command("train")
@click.option(
    "--embeddings-dir", "-i",
    default="./data/ml/embeddings",
    show_default=True,
    type=click.Path(exists=True),
)
@click.option(
    "--output-dir", "-o",
    default="./results/ml",
    show_default=True,
)
@click.option("--lr", default=1e-3, show_default=True, type=float)
@click.option("--dropout", default=0.3, show_default=True, type=float)
@click.option("--pos-weight", default=3.0, show_default=True, type=float)
@click.option("--batch-size", default=256, show_default=True)
@click.option("--max-epochs", default=50, show_default=True)
@click.option("--patience", default=7, show_default=True)
@click.option("--device", default="auto", show_default=True)
def ml_train(
    embeddings_dir: str,
    output_dir: str,
    lr: float,
    dropout: float,
    pos_weight: float,
    batch_size: int,
    max_epochs: int,
    patience: int,
    device: str,
):
    """Train the MLP classifier on pre-computed ESM-2 embeddings."""
    import numpy as np
    from paleoamp.ml.model import AMPClassifier, TrainConfig, train

    emb_dir = Path(embeddings_dir)
    X_train = np.load(emb_dir / "embeddings_train.npy")
    y_train = np.load(emb_dir / "labels_train.npy")
    X_val   = np.load(emb_dir / "embeddings_test.npy")
    y_val   = np.load(emb_dir / "labels_test.npy")

    _console.print(
        f"[bold]Training data:[/] {X_train.shape[0]:,} sequences  "
        f"(positives: {y_train.sum():,.0f})"
    )
    _console.print(
        f"[bold]Val data:[/]      {X_val.shape[0]:,} sequences  "
        f"(positives: {y_val.sum():,.0f})"
    )

    config = TrainConfig(
        input_dim=X_train.shape[1],
        dropout=dropout,
        lr=lr,
        pos_weight=pos_weight,
        batch_size=batch_size,
        max_epochs=max_epochs,
        patience=patience,
        device=device,
    )

    ckpt_dir = Path(output_dir) / "checkpoints"
    model, result = train(X_train, y_train, X_val, y_val, config, ckpt_dir)
    _console.print(
        f"\n[bold green]Training complete.[/]  "
        f"Best epoch: {result.best_epoch}  val_AUPR: {result.best_val_aupr:.4f}"
    )


@ml_group.command("score")
@click.option(
    "--orfs-dir", "-i",
    default="./results/orfs_hydrophobic",
    show_default=True,
    help="Directory of per-sample .faa files to score.",
)
@click.option(
    "--checkpoint", "-m",
    default="./results/ml/checkpoints/best_model.pt",
    show_default=True,
    type=click.Path(exists=True),
    help="Trained MLP checkpoint (.pt file).",
)
@click.option(
    "--output-dir", "-o",
    default="./results/ml/scores",
    show_default=True,
)
@click.option(
    "--threshold", default=0.5, show_default=True, type=float,
    help="Probability threshold to call a sequence AMP-positive.",
)
@click.option(
    "--latest-only/--all-samples", default=True, show_default=True,
    help="Score only the most recently modified sample (--latest-only) or all samples (--all-samples).",
)
@click.option("--batch-size", default=64, show_default=True)
@click.option("--device", default="auto", show_default=True)
def ml_score(
    orfs_dir: str,
    checkpoint: str,
    output_dir: str,
    threshold: float,
    latest_only: bool,
    batch_size: int,
    device: str,
):
    """
    Score ORF sequences with the trained AMP classifier.

    By default (--latest-only) only the most recently downloaded/assembled
    sample is scored so the output reflects the current pipeline run.
    Use --no-latest-only to score all samples in --orfs-dir.
    """
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    from paleoamp.ml.embed import embed_sequences
    from paleoamp.ml.model import load_model, predict_proba
    from rich.table import Table
    from rich import box

    orfs_path = Path(orfs_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if latest_only:
        faa_files = [_latest_sample_orfs(orfs_path)]
        _console.print(f"[bold]Scoring latest sample:[/] {faa_files[0].parent.name}")
    else:
        faa_files = [
            d / f"{d.name}.faa"
            for d in sorted(orfs_path.iterdir())
            if d.is_dir() and (d / f"{d.name}.faa").exists()
        ]

    dev = ("cuda" if __import__("torch").cuda.is_available() else "cpu") if device == "auto" else device
    model = load_model(Path(checkpoint), device=dev)

    table = Table(title="PaleoAMP — AMP Classifier Scores", box=box.ROUNDED, highlight=True)
    table.add_column("Sample", style="cyan")
    table.add_column("ORFs", justify="right")
    table.add_column(f"Score ≥ {threshold}", justify="right")
    table.add_column("Top hit", style="dim")

    for faa in faa_files:
        sample_id = faa.parent.name
        records = list(SeqIO.parse(str(faa), "fasta"))
        seqs = [str(r.seq).upper().rstrip("*") for r in records]
        ids  = [r.id for r in records]

        if not seqs:
            table.add_row(sample_id, "0", "0", "[dim]SKIPPED — no sequences[/]")
            continue

        embeddings = embed_sequences(seqs, batch_size=batch_size, device=dev)
        scores = predict_proba(model, embeddings, batch_size=batch_size, device=dev)

        df_out = pd.DataFrame({"seq_id": ids, "sequence": seqs, "amp_score": scores})
        df_out = df_out.sort_values("amp_score", ascending=False).reset_index(drop=True)

        tsv_path = out_path / f"{sample_id}_amp_scores.tsv"
        df_out.to_csv(tsv_path, sep="\t", index=False)

        n_pos = (scores >= threshold).sum()
        top = df_out.iloc[0]["seq_id"] if len(df_out) else "—"
        table.add_row(sample_id, f"{len(seqs):,}", f"{n_pos:,}", top)

    _console.print(table)


@ml_group.command("novelty")
@click.option(
    "--scores-tsv", "-i",
    default=None,
    help=(
        "Scores TSV from 'paleoamp ml score'. "
        "Defaults to the most recent file in --scores-dir."
    ),
)
@click.option(
    "--scores-dir",
    default="./results/ml/scores",
    show_default=True,
)
@click.option(
    "--amp-db", "-d",
    default="./data/amp_databases/merged_amp_db.fasta",
    show_default=True,
    type=click.Path(exists=True),
    help="AMP reference FASTA to screen against (from 'paleoamp db merge').",
)
@click.option(
    "--output-dir", "-o",
    default="./results/ml/novelty",
    show_default=True,
)
@click.option(
    "--score-threshold", default=0.5, show_default=True, type=float,
    help="Minimum AMP score to include a sequence as a candidate.",
)
@click.option(
    "--pident", default=40.0, show_default=True, type=float,
    help="Minimum % identity to consider a sequence known (not novel).",
)
@click.option(
    "--qcov", default=80.0, show_default=True, type=float,
    help="Minimum query coverage % to consider a sequence known.",
)
@click.option("--threads", "-t", default=4, show_default=True)
def ml_novelty_cmd(
    scores_tsv: str | None,
    scores_dir: str,
    amp_db: str,
    output_dir: str,
    score_threshold: float,
    pident: float,
    qcov: float,
    threads: int,
):
    """
    Screen high-scoring AMP candidates for novelty using MMseqs2.

    Candidates with no match to known AMPs at ≥ --pident identity and
    ≥ --qcov query coverage are reported as candidate ancient novel AMPs.
    """
    from paleoamp.ml.novelty import screen_novelty
    from rich.table import Table
    from rich import box

    # Resolve scores TSV — default to most recent file in scores-dir
    if scores_tsv:
        tsv_path = Path(scores_tsv)
    else:
        candidates = sorted(Path(scores_dir).glob("*_amp_scores.tsv"), key=lambda p: p.stat().st_mtime)
        if not candidates:
            _console.print(f"[red]No score files found in {scores_dir}[/]")
            sys.exit(1)
        tsv_path = candidates[-1]

    sample_id = tsv_path.stem.replace("_amp_scores", "")
    _console.print(f"[bold]Novelty screen:[/] {tsv_path.name}  (threshold ≥ {score_threshold})")

    result = screen_novelty(
        scores_tsv=tsv_path,
        amp_db_fasta=Path(amp_db),
        output_dir=Path(output_dir),
        sample_id=sample_id,
        score_threshold=score_threshold,
        pident_threshold=pident,
        qcov_threshold=qcov,
        threads=threads,
    )

    table = Table(title="PaleoAMP — Novelty Screen", box=box.ROUNDED, highlight=True)
    table.add_column("Sample", style="cyan")
    table.add_column("Candidates", justify="right")
    table.add_column("Known hits", justify="right")
    table.add_column("Novel", justify="right", style="bold green")
    table.add_row(
        result.sample_id,
        str(result.n_candidates),
        str(result.n_hits),
        str(result.n_novel),
    )
    _console.print(table)

    if result.n_novel:
        _console.print(f"\nTop novel candidates (by AMP score):")
        top = result.novel.head(10)
        for _, row in top.iterrows():
            _console.print(f"  [cyan]{row['seq_id']}[/]  score={row['amp_score']:.3f}  {row['sequence'][:40]}{'…' if len(row['sequence']) > 40 else ''}")
        _console.print(f"\nFull list → [cyan]{result.tsv_path}[/]")
    else:
        _console.print("[yellow]No novel candidates found above threshold.[/]")


@ml_group.command("validate")
@click.option(
    "--novel-tsv", "-i",
    default=None,
    help="Novel candidates TSV from 'paleoamp ml novelty'. Defaults to most recent.",
)
@click.option(
    "--novelty-dir",
    default="./results/ml/novelty",
    show_default=True,
)
@click.option(
    "--gff-dir",
    default="./results/orfs",
    show_default=True,
    help="Directory of per-sample Prodigal .gff files.",
)
@click.option(
    "--dataset", "-d",
    default="./data/ml/dataset_split.csv",
    show_default=True,
    type=click.Path(exists=True),
    help="Training dataset CSV used to fit the AAC second classifier.",
)
@click.option(
    "--output-dir", "-o",
    default="./results/ml/validated",
    show_default=True,
)
@click.option(
    "--aac-threshold", default=0.4, show_default=True, type=float,
    help="Minimum AAC classifier score to count as a positive vote.",
)
def ml_validate_cmd(
    novel_tsv: str | None,
    novelty_dir: str,
    gff_dir: str,
    dataset: str,
    output_dir: str,
    aac_threshold: float,
):
    """
    Validate novel AMP candidates through three independent layers:

    \b
    1. Physicochemical filters (charge, amphipathicity, length, complexity)
    2. Prodigal GFF partial-ORF flags
    3. AAC logistic regression (second classifier, independent of ESM-2)

    Verdict: PASS (3-4 criteria met), WARN (2), FAIL (0-1), DUPLICATE.
    """
    from paleoamp.ml.validate import validate_candidates
    from rich.table import Table
    from rich import box

    if novel_tsv:
        tsv_path = Path(novel_tsv)
    else:
        candidates = sorted(
            Path(novelty_dir).glob("*_novel.tsv"),
            key=lambda p: p.stat().st_mtime,
        )
        if not candidates:
            _console.print(f"[red]No novel TSV files found in {novelty_dir}[/]")
            sys.exit(1)
        tsv_path = candidates[-1]

    sample_id = tsv_path.stem.replace("_novel", "")
    gff_path = Path(gff_dir) / sample_id / f"{sample_id}.gff"
    if not gff_path.exists():
        _console.print(f"[red]GFF not found: {gff_path}[/]")
        sys.exit(1)

    _console.print(f"[bold]Validating:[/] {tsv_path.name}")

    df = validate_candidates(
        novel_tsv=tsv_path,
        gff_path=gff_path,
        dataset_csv=Path(dataset),
        output_dir=Path(output_dir),
        aac_threshold=aac_threshold,
    )

    table = Table(title="PaleoAMP — Candidate Validation", box=box.ROUNDED, highlight=True)
    table.add_column("ORF", style="cyan")
    table.add_column("Len", justify="right")
    table.add_column("Charge", justify="right")
    table.add_column("H-moment", justify="right")
    table.add_column("Partial", justify="center")
    table.add_column("AAC", justify="right")
    table.add_column("ML score", justify="right")
    table.add_column("Verdict", justify="center")

    verdict_style = {
        "PASS": "bold green", "WARN": "yellow",
        "FAIL": "red", "DUPLICATE": "dim",
    }
    for _, row in df.iterrows():
        style = verdict_style.get(row["verdict"], "")
        table.add_row(
            row["seq_id"],
            str(int(row["length"])),
            f"{row['net_charge']:+.1f}",
            f"{row['hydrophobic_moment']:.3f}",
            row["partial_flag"],
            f"{row['aac_score']:.2f}",
            f"{row['amp_score']:.3f}",
            f"[{style}]{row['verdict']}[/]",
        )

    _console.print(table)

    counts = df["verdict"].value_counts()
    _console.print(
        f"\n[bold green]PASS: {counts.get('PASS', 0)}[/]  "
        f"[yellow]WARN: {counts.get('WARN', 0)}[/]  "
        f"[red]FAIL: {counts.get('FAIL', 0)}[/]  "
        f"[dim]DUPLICATE: {counts.get('DUPLICATE', 0)}[/]"
    )
    out = Path(output_dir) / f"{tsv_path.stem}_validated.tsv"
    _console.print(f"Full results → [cyan]{out}[/]")


# ---------------------------------------------------------------------------
# paleoamp qc
# ---------------------------------------------------------------------------

@main.group("qc")
def qc_group():
    """Assess ancient DNA quality of sequencing reads."""


@qc_group.command("assess")
@click.argument("input_path", metavar="INPUT")
@click.option(
    "--config", "-c",
    default=None,
    type=click.Path(exists=True),
    help="YAML config file (defaults to config/defaults.yaml).",
)
@click.option(
    "--output-dir", "-o",
    default="./results/qc",
    show_default=True,
    help="Directory to write JSON and TSV reports.",
)
@click.option(
    "--plot/--no-plot",
    default=True,
    show_default=True,
    help="Print ASCII damage plots for each file.",
)
@click.option(
    "--max-reads", "-n",
    default=None,
    type=int,
    help="Limit reads processed per file (useful for large files).",
)
@click.option(
    "--min-ct-rate",
    default=None,
    type=float,
    help="Override minimum 5′ C→T rate threshold.",
)
@click.option(
    "--min-ga-rate",
    default=None,
    type=float,
    help="Override minimum 3′ G→A rate threshold.",
)
def qc_assess(
    input_path: str,
    config: str | None,
    output_dir: str,
    plot: bool,
    max_reads: int | None,
    min_ct_rate: float | None,
    min_ga_rate: float | None,
):
    """
    Assess ancient DNA quality of FASTQ files.

    INPUT may be a single FASTQ/FASTQ.gz file or a directory containing them.
    """
    from paleoamp.qc.adna import QCThresholds, assess_fastq, assess_fastq_directory
    from paleoamp.qc.report import (
        print_summary_table,
        print_damage_plot,
        save_json_report,
        save_tsv_report,
    )

    cfg = _load_config(Path(config) if config else None)
    thresholds = QCThresholds.from_config(cfg)

    if min_ct_rate is not None:
        thresholds.min_ct_rate_5prime = min_ct_rate
    if min_ga_rate is not None:
        thresholds.min_ga_rate_3prime = min_ga_rate

    input_p = Path(input_path)
    if input_p.is_dir():
        _console.print(f"Assessing FASTQ files in [cyan]{input_p}[/] ...")
        profiles = assess_fastq_directory(input_p, thresholds=thresholds, max_reads=max_reads)
    elif input_p.is_file():
        _console.print(f"Assessing [cyan]{input_p}[/] ...")
        profiles = [assess_fastq(input_p, thresholds=thresholds, max_reads=max_reads)]
    else:
        _console.print(f"[red]Path not found: {input_p}[/]")
        sys.exit(1)

    print_summary_table(profiles)

    if plot:
        for p in profiles:
            print_damage_plot(p)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_json_report(profiles, out_dir / "qc_report.json")
    save_tsv_report(profiles, out_dir / "qc_report.tsv")

    n_pass = sum(1 for p in profiles if p.pass_verdict)
    _console.print(
        f"\n[bold]{n_pass}/{len(profiles)} file(s) passed quality thresholds.[/]"
    )
    sys.exit(0 if n_pass > 0 else 1)
