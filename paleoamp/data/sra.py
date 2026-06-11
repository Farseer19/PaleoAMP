"""Query and fetch ancient metagenomic datasets from NCBI SRA."""

from __future__ import annotations

import io
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests
from pysradb.sraweb import SRAweb  # still used by get_run_metadata

# NCBI E-utilities endpoints
_ESEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_EFETCH_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Safe eFetch batch size — above ~300 NCBI may drop the connection
_FETCH_BATCH = 200

# Standard columns for NCBI SRA runinfo CSV (no header row in eFetch response)
_RUNINFO_COLS = [
    "Run", "ReleaseDate", "LoadDate", "spots", "bases",
    "spots_with_mates", "avgLength", "size_MB", "AssemblyName",
    "download_path", "Experiment", "LibraryName", "LibraryStrategy",
    "LibrarySelection", "LibrarySource", "LibraryLayout",
    "InsertSize", "InsertDev", "Platform", "Model",
    "SRAStudy", "BioProject", "Study_Pubmed_id", "ProjectID",
    "Sample", "BioSample", "SampleType", "TaxID", "ScientificName",
    "SampleName", "g1k_pop_code", "source", "g1k_analysis_group",
    "Subject_ID", "Sex", "Disease", "Tumor", "Affection_Status",
    "Analyte_Type", "Histological_Type", "Body_Site",
    "CenterName", "Submission", "dbgap_study_accession",
    "Consent", "RunHash", "ReadHash",
]


def search_ancient_metagenomes(
    query: str,
    max_results: int = 200,
) -> pd.DataFrame:
    """
    Search SRA for datasets matching *query* and return a DataFrame.

    Uses NCBI E-utilities directly with paginated eFetch (runinfo CSV format)
    so that large result sets — e.g. --max-results 1000 — never trigger the
    ChunkedEncodingError that pysradb's single-shot eFetch causes.
    """
    try:
        return _batched_sra_search(query, max_results)
    except Exception as exc:
        raise RuntimeError(f"SRA search failed: {exc}") from exc


# ---------------------------------------------------------------------------
# Internal: paginated search
# ---------------------------------------------------------------------------

def _batched_sra_search(query: str, max_results: int) -> pd.DataFrame:
    # Step 1 — eSearch: get UIDs and register with NCBI history server so
    # subsequent eFetch calls can page through the full result set.
    r = requests.get(_ESEARCH_URL, params={
        "db": "sra",
        "term": query,
        "retmax": min(max_results, 100_000),
        "retmode": "json",
        "usehistory": "y",
    }, timeout=30)
    r.raise_for_status()

    esearch = r.json()["esearchresult"]
    total = int(esearch["count"])
    if total == 0:
        return pd.DataFrame()

    web_env  = esearch["webenv"]
    query_key = esearch["querykey"]
    to_fetch = min(max_results, total)

    print(
        f"Found {total:,} matching SRA records; fetching up to {to_fetch:,}.",
        file=sys.stderr,
    )

    # Step 2 — eFetch in batches using runinfo CSV (lightweight vs full XML).
    frames: list[pd.DataFrame] = []
    for start in range(0, to_fetch, _FETCH_BATCH):
        size = min(_FETCH_BATCH, to_fetch - start)
        r = requests.get(_EFETCH_URL, params={
            "db": "sra",
            "query_key": query_key,
            "WebEnv": web_env,
            "retstart": start,
            "retmax": size,
            "rettype": "runinfo",
            "retmode": "csv",
        }, timeout=120)
        r.raise_for_status()
        text = r.text.strip()
        if not text:
            continue
        try:
            frames.append(_read_runinfo_csv(text))
        except Exception:
            continue

    if not frames:
        return pd.DataFrame()

    df = pd.concat(frames, ignore_index=True)
    if "Run" in df.columns:
        df = df.drop_duplicates(subset=["Run"])

    # Step 3 — Enrich with study titles (not in runinfo CSV).
    study_accs = df["SRAStudy"].dropna().unique().tolist()
    titles = _fetch_study_titles(study_accs)
    df["study_title"] = df["SRAStudy"].map(titles)

    return _normalise_runinfo(df)


def _read_runinfo_csv(text: str) -> pd.DataFrame:
    """
    Parse an NCBI SRA runinfo CSV response.

    NCBI's eFetch does not include a header row, so the standard column
    names are supplied explicitly. Any extra/missing columns are handled
    gracefully by truncating or padding the name list to match.
    """
    # Count columns from first line to handle NCBI format variations
    first_line = text.splitlines()[0]
    n_cols = len(first_line.split(","))
    col_names = _RUNINFO_COLS[:n_cols]
    # Pad with generic names if NCBI adds extra columns
    if n_cols > len(col_names):
        col_names += [f"extra_{i}" for i in range(n_cols - len(col_names))]
    return pd.read_csv(io.StringIO(text), names=col_names, header=None)


def _fetch_study_titles(study_accessions: list[str]) -> dict[str, str]:
    """
    Return a dict of SRP/ERP accession → study title using NCBI eSummary.

    Batches 20 studies per request to avoid per-query overhead while keeping
    response sizes manageable.
    """
    if not study_accessions:
        return {}

    titles: dict[str, str] = {}
    batch_size = 20

    for i in range(0, len(study_accessions), batch_size):
        batch = study_accessions[i : i + batch_size]
        term = " OR ".join(f'"{acc}"[Accession]' for acc in batch)
        try:
            # eSearch to convert accessions → UIDs
            r = requests.get(_ESEARCH_URL, params={
                "db": "sra",
                "term": term,
                "retmax": len(batch) * 10,
                "retmode": "json",
                "usehistory": "y",
            }, timeout=20)
            r.raise_for_status()
            esearch = r.json()["esearchresult"]
            if int(esearch["count"]) == 0:
                continue

            # eSummary — fetch one representative record per UID, parse expxml
            r = requests.get(_ESUMMARY_URL, params={
                "db": "sra",
                "query_key": esearch["querykey"],
                "WebEnv": esearch["webenv"],
                "retmax": len(batch) * 10,
                "retmode": "json",
            }, timeout=20)
            r.raise_for_status()

            for uid, data in r.json().get("result", {}).items():
                if uid == "uids":
                    continue
                expxml = data.get("expxml", "")
                if not expxml:
                    continue
                try:
                    root = ET.fromstring(f"<root>{expxml}</root>")
                    study_el = root.find(".//Study")
                    if study_el is not None:
                        acc = study_el.get("acc", "")
                        name = study_el.get("name", "")
                        if acc and name and acc not in titles:
                            titles[acc] = name
                except ET.ParseError:
                    continue
        except Exception:
            continue

    return titles


def _normalise_runinfo(df: pd.DataFrame) -> pd.DataFrame:
    """Rename runinfo CSV columns to our standard schema and add ENA URLs."""
    rename = {
        "Run":            "run_accession",
        "SRAStudy":       "study_accession",
        "Experiment":     "experiment_accession",
        "LibraryStrategy":"library_strategy",
        "LibrarySource":  "library_source",
        "LibraryLayout":  "library_layout",
        "Platform":       "instrument",
        "Model":          "instrument_model",
        "spots":          "total_spots",
        "bases":          "total_bases",
        "size_MB":        "total_size_mb",
        "TaxID":          "organism_taxid",
        "ScientificName": "organism_name",
        "SampleName":     "sample_title",
        "BioSample":      "biosample",
        "BioProject":     "bioproject",
        "CenterName":     "center_name",
        "LibraryName":    "library_name",
        "collection_date":"collection_date",
        "geo_loc_name_country": "geo_loc_name",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    if "run_accession" in df.columns:
        df["ena_fastq_http"] = df["run_accession"].apply(
            lambda acc: _ena_ftp_url(str(acc))
            if pd.notna(acc) and str(acc)[:3] in ("ERR", "SRR", "DRR")
            else ""
        )

    return df


def filter_for_ancient_indicators(df: pd.DataFrame) -> pd.DataFrame:
    """
    Heuristic filter: keep rows whose metadata fields suggest ancient DNA.

    Searches across title, description, biome, and material columns.
    """
    if df.empty:
        return df

    ancient_keywords = [
        "ancient", "paleo", "fossil", "permafrost", "archaeological",
        "holocene", "pleistocene", "ice core", "mummy", "archaic",
        "subfossil", "historic", "medieval", "burial", "extinct",
        "coprolite", "environmental dna", "late quaternary",
    ]
    pattern = "|".join(ancient_keywords)

    search_cols = [
        c for c in [
            "study_title", "experiment_title", "experiment_desc",
            "sample_title", "organism_name", "library_name",
            "environment (biome)", "environment (material)",
            "geo_loc_name", "isolation_source",
        ]
        if c in df.columns
    ]

    if not search_cols:
        return df

    combined = df[search_cols].fillna("").astype(str).agg(" ".join, axis=1)
    mask = combined.str.lower().str.contains(pattern, na=False)
    return df[mask].reset_index(drop=True)


def filter_for_assembly_size(
    df: pd.DataFrame,
    min_spots: int | None = None,
    min_bases: int | None = None,
    min_size_mb: float | None = None,
) -> pd.DataFrame:
    """
    Keep only runs with enough data for a meaningful metagenome assembly.

    Any combination of thresholds may be supplied; a row must satisfy ALL
    provided thresholds to be retained.  Rows where the relevant column is
    missing or non-numeric are dropped when a threshold is active.

    Typical guidance for aDNA metagenomes (short reads, ~40-80 bp):
      - min_spots  ≥ 500_000  — ~500 K reads, low-end workable assembly
      - min_bases  ≥ 5×10^8   — 500 Mbp raw bases
      - min_size_mb ≥ 500     — ~500 MB compressed FASTQ
    """
    if df.empty or (min_spots is None and min_bases is None and min_size_mb is None):
        return df

    mask = pd.Series(True, index=df.index)

    def _numeric(col: str) -> pd.Series:
        return pd.to_numeric(df[col], errors="coerce") if col in df.columns else None

    if min_spots is not None:
        s = _numeric("total_spots")
        if s is not None:
            mask &= s >= min_spots
        else:
            mask &= False

    if min_bases is not None:
        b = _numeric("total_bases")
        if b is not None:
            mask &= b >= min_bases
        else:
            mask &= False

    if min_size_mb is not None:
        m = _numeric("total_size_mb")
        if m is not None:
            mask &= m >= min_size_mb
        else:
            mask &= False

    return df[mask].reset_index(drop=True)


def filter_for_microbial(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep rows that are likely to be primarily microbial metagenomes.

    Applies two criteria in order:
    1. library_source must be METAGENOMIC (drops amplicon, WGS of isolates, etc.)
    2. At least one positive microbial keyword must appear in the metadata,
       AND the row must not match exclusion patterns for plant/animal eDNA studies.
    """
    if df.empty:
        return df

    # Criterion 1: METAGENOMIC library source
    if "library_source" in df.columns:
        df = df[df["library_source"].str.upper() == "METAGENOMIC"].copy()

    if df.empty:
        return df

    # Criterion 2a: positive microbial signal
    microbial_keywords = [
        "microbiome", "microbiota", "microbial", "bacteria", "archaea",
        "prokaryot", "dental calculus", "oral microbiome", "gut microbiome",
        "soil metagenome", "sediment metagenome", "permafrost metagenome",
        "marine metagenome", "freshwater metagenome", "hot spring",
        "gut", "oral", "skin microbiome", "rumen", "rhizosphere",
        "pathogen", "antimicrobial", "antibiotic resistance",
    ]

    # Criterion 2b: exclusion — studies explicitly targeting plants or animals
    # (environmental metagenomes that are clearly not microbial-focused)
    exclude_keywords = [
        "plant and animal", "plant edna", "animal edna",
        "vertebrate edna", "eukaryote", "herbarium",
        "insect metagenome", "zooplankton metagenome",
    ]

    search_cols = [
        c for c in [
            "study_title", "experiment_title", "experiment_desc",
            "sample_title", "organism_name", "library_name",
            "environment (biome)", "environment (material)",
            "isolation_source",
        ]
        if c in df.columns
    ]

    if not search_cols:
        return df

    combined = (
        df[search_cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower()
    )

    microbial_mask = combined.str.contains(
        "|".join(microbial_keywords), na=False, regex=True
    )
    exclude_mask = combined.str.contains(
        "|".join(exclude_keywords), na=False, regex=True
    )

    return df[microbial_mask & ~exclude_mask].reset_index(drop=True)


def download_via_ena(
    accessions: list[str],
    output_dir: Path,
    metadata_df: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """
    Download FASTQ.gz files directly from ENA HTTP URLs.

    Preferred over fasterq-dump when ENA URLs are available — no SRA Toolkit
    required.  Falls back to constructing the canonical ENA FTP path if the
    URL is not already in *metadata_df*.

    Returns accession → downloaded file path for successful downloads.
    """
    import requests

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Path] = {}

    # Build a lookup of accession → ENA HTTP URL from metadata if provided
    url_map: dict[str, list[str]] = {}
    if metadata_df is not None and not metadata_df.empty:
        for _, row in metadata_df.iterrows():
            acc = row.get("run_accession", "")
            if not acc:
                continue
            urls = []
            for col in ("ena_fastq_http", "ena_fastq_http_1", "ena_fastq_http_2"):
                val = row.get(col)
                val_str = str(val) if val is not None else ""
                if val_str.startswith("http"):
                    urls.append(val_str)
            if urls:
                url_map[acc] = urls

    for acc in accessions:
        urls = url_map.get(acc) or _ena_ftp_urls(acc)
        acc_dir = output_dir / acc
        acc_dir.mkdir(exist_ok=True)
        downloaded = []
        for url in urls:
            filename = url.split("/")[-1]
            dest = acc_dir / filename
            if dest.exists():
                print(f"[skip] {filename} already exists")
                downloaded.append(dest)
                continue
            print(f"[download] {acc}: {url}")
            try:
                with requests.get(url, stream=True, timeout=120) as resp:
                    if resp.status_code == 404:
                        continue  # silently skip; try next URL variant
                    resp.raise_for_status()
                    with open(dest, "wb") as fh:
                        for chunk in resp.iter_content(chunk_size=65536):
                            fh.write(chunk)
                size_mb = dest.stat().st_size / 1e6
                print(f"[ok] {filename} ({size_mb:.1f} MB)")
                downloaded.append(dest)
            except Exception as exc:
                print(f"[error] {acc}: {exc}", file=sys.stderr)
        if downloaded:
            results[acc] = acc_dir
    return results


def _ena_ftp_url(accession: str) -> str:
    """Construct the canonical ENA HTTP path for a run accession (single-end)."""
    return _ena_ftp_urls(accession)[0]


def _ena_ftp_urls(accession: str) -> list[str]:
    """
    Return candidate ENA HTTP URLs for *accession*, paired before single-end.

    ENA FTP path rules by accession length:
      9 chars (6-digit):  vol1/fastq/[pre6]/[acc]/
      10 chars (7-digit): vol1/fastq/[pre6]/[last3]/[acc]/
      11 chars (8-digit): vol1/fastq/[pre6]/0[last2]/[acc]/   (rare)
    """
    prefix = accession[:6]
    n = len(accession)
    if n == 9:
        base = f"vol1/fastq/{prefix}/{accession}"
    elif n == 10:
        base = f"vol1/fastq/{prefix}/{accession[-3:]}/{accession}"
    elif n == 11:
        base = f"vol1/fastq/{prefix}/0{accession[-2:]}/{accession}"
    else:
        base = f"vol1/fastq/{prefix}/{accession}"

    root = f"http://ftp.sra.ebi.ac.uk/{base}"
    return [
        f"{root}/{accession}_1.fastq.gz",
        f"{root}/{accession}_2.fastq.gz",
        f"{root}/{accession}.fastq.gz",
    ]


def download_runs(
    accessions: list[str],
    output_dir: Path,
    threads: int = 4,
    use_fasterq: bool = True,
) -> dict[str, Path]:
    """
    Download FASTQ files for *accessions* using fasterq-dump (preferred) or
    fastq-dump as fallback.

    Returns a mapping of accession → output directory path for successfully
    downloaded runs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tool = "fasterq-dump" if use_fasterq else "fastq-dump"
    _check_sra_tools(tool)

    results: dict[str, Path] = {}
    for acc in accessions:
        acc_dir = output_dir / acc
        acc_dir.mkdir(exist_ok=True)
        cmd = [
            tool, acc,
            "--outdir", str(acc_dir),
            "--threads", str(threads),
            "--progress",
        ]
        if tool == "fasterq-dump":
            cmd += ["--split-3"]
        try:
            subprocess.run(cmd, check=True)
            results[acc] = acc_dir
        except subprocess.CalledProcessError as exc:
            print(f"[warn] Failed to download {acc}: {exc}", file=sys.stderr)
    return results


def _check_sra_tools(tool: str) -> None:
    if subprocess.run(["which", tool], capture_output=True).returncode != 0:
        raise EnvironmentError(
            f"{tool} not found. Install SRA Toolkit: "
            "https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
        )


def get_run_metadata(accessions: list[str]) -> pd.DataFrame:
    """Fetch run-level metadata for a list of SRA accessions."""
    db = SRAweb()
    frames = []
    for acc in accessions:
        try:
            meta = db.sra_metadata(acc, detailed=True)
            if meta is not None and not meta.empty:
                frames.append(meta)
        except Exception as exc:
            print(f"[warn] Metadata fetch failed for {acc}: {exc}", file=sys.stderr)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
