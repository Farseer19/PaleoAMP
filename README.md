# PaleoAMP

A command-line pipeline for discovering novel antimicrobial peptide (AMP) genes
in ancient and extinct microbial samples from paleogenomic sequencing data.

Ancient metagenomes preserve DNA from microbes no longer captured in modern
culture collections or databases. PaleoAMP mines that dark matter by running
raw reads through assembly, ORF prediction, and an ESM-2–powered classifier to
surface candidate AMPs that have no close relative in any known AMP database.

---

## Overview

```
Raw FASTQ reads (SRA)
        │
        ▼
  aDNA Quality Control          ← paleoamp qc assess
        │  (filter degraded / low-quality samples)
        ▼
  Metagenomic Assembly           ← paleoamp assemble run   (MEGAHIT)
        │  (contig set per sample)
        ▼
  ORF Prediction                 ← paleoamp predict orfs   (Prodigal)
        │  (protein FASTA + GFF per sample)
        ▼
  Hydrophobicity Pre-filter      ← paleoamp predict hydrophobic
        │  (retain AMP-like physicochemical profile)
        ▼
  ESM-2 AMP Classifier           ← paleoamp ml score
        │  (probability score for each ORF)
        ▼
  Novelty Screen (MMseqs2)       ← paleoamp ml novelty
        │  (remove sequences matching known AMPs)
        ▼
  Multi-criteria Validation      ← paleoamp ml validate
        │  (charge, amphipathicity, partial-ORF, AAC classifier)
        ▼
  Novel AMP Candidates TSV
```

The ML classifier is trained once on curated positive (APD, DRAMP general,
DRAMP specific, UniProt reviewed AMPs) and negative (UniProt non-AMP)
sequences using frozen ESM-2 embeddings.
Training commands are described in the [Classifier Training](#classifier-training) section.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Farseer19/PaleoAMP.git
cd PaleoAMP
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate paleoamp
```

This installs all Python packages and the required bioinformatics tools
(MEGAHIT, Prodigal, MMseqs2) from the bioconda channel.

### 3. Install the PaleoAMP package

```bash
pip install -e .
```

Verify the install:

```bash
paleoamp --version
```

---

## Typical Command Flow

### Step 1 — Set up the AMP reference database

Download the three source databases and merge them into a single
non-redundant FASTA. This file is used by `paleoamp ml novelty` as the
MMseqs2 reference to determine whether a candidate AMP is already known.

```bash
paleoamp db download
paleoamp db merge
```

Databases are saved to `data/amp_databases/` by default. Pass `--force` to
re-download if files already exist.

---

### Step 2 — Find and download ancient metagenome reads

#### `paleoamp sra search`

Queries NCBI SRA via the E-utilities API and returns a table of matching
runs. Results are printed to the terminal and optionally saved to a TSV file.

```bash
paleoamp sra search "permafrost metagenome ancient" \
    --microbial-only \
    --min-reads 500000 \
    --output candidates.tsv
```

**Filters:**

`--ancient-filter` *(on by default, disable with `--no-ancient-filter`)*
Scans every metadata field (study title, experiment description, sample
name, organism, geographic location, isolation source, and others) for
keywords that indicate an ancient or paleogenomic context:
`ancient`, `paleo`, `fossil`, `permafrost`, `archaeological`, `holocene`,
`pleistocene`, `ice core`, `mummy`, `archaic`, `subfossil`, `historic`,
`medieval`, `burial`, `extinct`, `coprolite`, `environmental dna`,
`late quaternary`.
This is a heuristic — it catches clearly ancient datasets but may miss
studies that don't use these exact terms in their metadata. Disable it
with `--no-ancient-filter` if you want every result back and intend to
inspect them manually.

`--microbial-only`
Applies two stricter criteria on top of the ancient filter:
1. The run's `library_source` field must be `METAGENOMIC` — this removes
   whole-genome sequencing of isolates, amplicon surveys (16S), and
   eukaryotic WGS.
2. At least one microbial context keyword must appear in the metadata
   (`microbiome`, `bacteria`, `archaea`, `dental calculus`, `soil
   metagenome`, `permafrost metagenome`, etc.) and no plant/animal eDNA
   exclusion keywords (`plant eDNA`, `vertebrate eDNA`, `herbarium`,
   `insect metagenome`, etc.). This keeps the focus on prokaryotic
   communities where novel AMPs are most likely to be found.

`--min-reads N`
Drops runs with fewer than N total read pairs (`total_spots`). Ancient
metagenomes typically have low microbial read counts after host and
environmental DNA subtraction — ≥ 500,000 reads is a practical lower
bound for a workable assembly.

`--min-bases BP`
Drops runs with fewer than BP total sequenced bases. An alternative to
`--min-reads` when comparing libraries with different read lengths.

`--min-size-mb MB`
Drops runs whose reported compressed file size is below MB megabytes.
Useful as a quick proxy when spot/base counts are not populated in the
metadata.

`--max-results N` *(default: 50)*
Maximum number of SRA records to return. Increase for broad surveys.

---

#### `paleoamp sra merge-candidates`

After running multiple searches with different query terms, merge the
resulting TSV files into one deduplicated table.

```bash
paleoamp sra merge-candidates search1.tsv search2.tsv -o all_candidates.tsv
```

Deduplicates on `run_accession` (keeping the first occurrence) and adds
a `found_in` column showing which source files each run appeared in.
Useful for knowing which accessions turned up across several independent
searches, which can be a signal of relevance.

---

#### `paleoamp sra download`

Downloads FASTQ files for one or more accessions.

```bash
paleoamp sra download ERR6458500 SRR33371653 SRR35641057
```

`--via-ena` *(default)*
Downloads directly from the European Nucleotide Archive over HTTP — no
SRA Toolkit required. Constructs the canonical ENA FTP path from the
accession format and tries paired-end files first (`_1.fastq.gz`,
`_2.fastq.gz`), falling back to single-end (`accession.fastq.gz`).

`--via-sratools`
Uses `fasterq-dump` from the NCBI SRA Toolkit instead. Required for
accessions not mirrored on ENA (rare, but some very new submissions).

Files are saved to `data/reads/<accession>/`.

---

### Step 3 — Assess ancient DNA quality

PaleoAMP checks for the characteristic cytosine deamination signature of
ancient DNA (elevated C→T at the 5′ end, G→A at the 3′ end). Only samples
that pass both damage thresholds are forwarded to assembly.

```bash
paleoamp qc assess data/reads/ --output-dir results/qc/
```

Pass a single file or a directory. The QC report (`results/qc/qc_report.json`)
is consumed automatically by `paleoamp assemble run`.

---

### Step 4 — Assemble reads into contigs

Runs MEGAHIT on each sample that passed aDNA QC. Samples that failed QC are
automatically skipped.

```bash
paleoamp assemble run \
    --reads-dir data/reads/ \
    --qc-report results/qc/qc_report.json \
    --output-dir results/assembly/ \
    --threads 8
```

---

### Step 5 — Predict open reading frames

Runs Prodigal in metagenome mode on each assembled sample. Outputs a protein
FASTA (`.faa`) and gene coordinate file (`.gff`) per sample.

```bash
paleoamp predict orfs \
    --assembly-dir results/assembly/ \
    --output-dir results/orfs/
```

---

### Step 6 — Hydrophobicity pre-filter

Retains only ORFs with an AMP-like physicochemical profile (Kyte-Doolittle
score and hydrophobic residue fraction) before the expensive ML step.

```bash
paleoamp predict hydrophobic \
    --orfs-dir results/orfs/ \
    --output-dir results/orfs_hydrophobic/
```

Default thresholds: mean KD ≥ 0.0, hydrophobic fraction ≥ 30 %.

---

### Step 7 — Score ORFs with the AMP classifier

Runs the trained ESM-2 + MLP model on the filtered ORFs and writes a
probability score for each sequence.

```bash
paleoamp ml score \
    --orfs-dir results/orfs_hydrophobic/ \
    --checkpoint results/ml/checkpoints/best_model.pt \
    --output-dir results/ml/scores/
```

Use `--all-samples` to score every sample in `--orfs-dir` in one pass
instead of only the most recently modified one.

---

### Step 8 — Screen for novelty

Compares high-scoring candidates against the merged AMP database using
MMseqs2. Sequences with no hit at ≥ 40 % identity and ≥ 80 % query coverage
are reported as candidate novel ancient AMPs.

```bash
paleoamp ml novelty \
    --amp-db data/amp_databases/merged_amp_db.fasta \
    --output-dir results/ml/novelty/
```

---

### Step 9 — Validate candidates

Three independent validation layers reduce false positives before final
reporting:

1. **Physicochemical filters** — net charge, hydrophobic moment, length,
   low-complexity masking, duplicate removal
2. **Prodigal GFF partial-ORF check** — flags incomplete gene calls at contig
   edges
3. **AAC logistic regression** — a second classifier trained on amino acid
   composition, fully independent of the ESM-2 model

```bash
paleoamp ml validate \
    --gff-dir results/orfs/ \
    --dataset data/ml/dataset_split.csv \
    --output-dir results/ml/validated/
```

Each candidate receives a verdict: **PASS** (3–4 criteria met), **WARN** (2),
**FAIL** (0–1), or **DUPLICATE**.

---

## Classifier Training

Run these commands once to train the ESM-2 + MLP classifier. Skip this if
you are using the pre-trained checkpoint at `results/ml/checkpoints/best_model.pt`.

```bash
# 1. Download training sequences and build the labelled dataset
paleoamp ml collect --output-dir data/ml/

# 2. Cluster at 40 % identity and split into train/test without data leakage
paleoamp ml cluster --dataset data/ml/dataset.csv

# 3. Generate frozen ESM-2 embeddings (GPU recommended; ~15 min on CPU)
paleoamp ml embed --dataset data/ml/dataset_split.csv

# 4. Train the MLP (typically < 50 epochs with early stopping)
paleoamp ml train --embeddings-dir data/ml/embeddings/ --output-dir results/ml/
```

---

## Configuration

Default thresholds live in `config/defaults.yaml`. Key settings:

| Section | Parameter | Default | Description |
|---|---|---|---|
| `quality` | `min_mean_phred` | 20 | Minimum mean base quality |
| `quality` | `min_read_length` | 30 bp | Minimum read length |
| `adna_damage` | `min_ct_rate_5prime` | 0.05 | 5′ C→T damage threshold |
| `adna_damage` | `min_ga_rate_3prime` | 0.05 | 3′ G→A damage threshold |
| `assembly` | `min_contig_len` | 200 bp | Minimum assembled contig length |

Override thresholds at runtime with flags such as `--min-ct-rate` and
`--min-ga-rate` on `paleoamp qc assess`, or edit `config/defaults.yaml`
to change the project-wide defaults.

---

## Output Files

| Path | Contents |
|---|---|
| `results/qc/qc_report.json` | Per-sample aDNA QC verdicts (used by assembler) |
| `results/qc/qc_report.tsv` | Same report in tabular form |
| `results/assembly/<id>/final.contigs.fa` | Assembled contigs per sample |
| `results/orfs/<id>/<id>.faa` | Predicted proteins per sample |
| `results/orfs/<id>/<id>.gff` | Gene coordinates per sample |
| `results/orfs_hydrophobic/<id>/<id>.faa` | Hydrophobicity-filtered proteins |
| `results/ml/scores/<id>_amp_scores.tsv` | AMP probability score per ORF |
| `results/ml/novelty/<id>_novel.tsv` | Candidate novel AMPs |
| `results/ml/novelty/<id>_known_hits.tsv` | ORFs matched to known AMPs |
| `results/ml/validated/<id>_novel_validated.tsv` | Final validated candidates |
