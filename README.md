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
        │  (C→T / G→A damage thresholds)
        ▼
  Adapter Trimming (optional)   ← paleoamp qc trim   (fastp)
        │
        ▼
  Metagenomic Assembly           ← paleoamp assemble run   (MEGAHIT)
        │  (contig set per sample)
        ▼
  ORF Prediction                 ← paleoamp predict orfs   (Prodigal)
        │  (protein FASTA + GFF per sample)
        ▼
  Amphipathicity Pre-filter      ← paleoamp predict hydrophobic
        │  (Eisenberg hydrophobic moment µH ≥ 0.10)
        ▼
  ESM-2 AMP Classifier           ← paleoamp ml score
        │  (probability score per ORF)
        ▼
  Novelty Screen (MMseqs2)       ← paleoamp ml novelty
        │  (remove sequences matching known AMPs at ≥ 40% identity)
        ▼
  Multi-criteria Validation      ← paleoamp ml validate
        │  (physicochemical, Boman index, partial-ORF, AAC classifier)
        ▼
  Novel AMP Candidates TSV
```

The ML classifier is trained on positive sequences (APD, DRAMP, UniProt AMPs,
including bacterial and archaeal subsets) and negative sequences (UniProt
non-AMP, transmembrane segments, signal peptides, shuffled positives).
Training commands are described in [Classifier Training](#classifier-training).

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
(MEGAHIT, Prodigal, MMseqs2, fastp) from the bioconda channel.

### 3. Install the PaleoAMP package

```bash
pip install -e .
```

---

## Typical Command Flow

### Step 1 — Set up the AMP reference database

Download AMP source databases to `data/amp_sources/` (shared with the ML
training pipeline) and merge them into a single non-redundant FASTA used by
the novelty screen.

```bash
paleoamp db download
paleoamp db merge
```

Downloads APD2024a, DRAMP 3.0 general and specific, UniProt reviewed AMPs
(KW-0929 including bacterial taxonomy:2 and archaeal taxonomy:2157 subsets),
dbAMP v2.0, DBAASP, and BACTIBASE. Sequences containing non-standard amino
acids (B, O, U, Z, alignment hyphens) are automatically excluded from the
merged file to prevent silent MMseqs2 failures.

---

### Step 2 — Find and download ancient metagenome reads

#### `paleoamp sra search`

Queries NCBI SRA and returns matching runs.

```bash
paleoamp sra search "permafrost metagenome ancient" \
    --microbial-only \
    --min-reads 500000 \
    --output candidates.tsv
```

**Filters:**

`--ancient-filter` *(on by default)*
Scans metadata fields for keywords indicating an ancient context:
`ancient`, `paleo`, `fossil`, `permafrost`, `archaeological`, `holocene`,
`pleistocene`, `ice core`, `mummy`, `archaic`, `subfossil`, `historic`,
`medieval`, `burial`, `extinct`, `coprolite`, `environmental dna`,
`late quaternary`.

`--microbial-only`
Requires `library_source = METAGENOMIC` and at least one microbial context
keyword. Excludes plant/animal eDNA studies.

`--min-reads N`
Drops runs with fewer than N total reads. ≥ 500,000 is a practical floor
for a workable aDNA assembly.

`--max-results N` *(default: 50)*

---

#### `paleoamp sra download`

```bash
paleoamp sra download ERR6458500 SRR33371653 SRR35641057
```

Downloads directly from ENA over HTTP by default (`--via-ena`). Files are
saved to `data/reads/<accession>/`.

---

### Step 3 — Assess ancient DNA quality

Checks for the characteristic cytosine deamination signature (C→T at 5′,
G→A at 3′). Generates a bar chart of damage rates per terminal position.

```bash
paleoamp qc assess data/reads/ --output-dir results/qc/
```

---

### Step 4 — Adapter trimming (optional)

Many ancient eDNA libraries deposited in SRA are pre-trimmed. Run this step
if you observe adapter contamination in the QC report.

```bash
paleoamp qc trim data/reads/
```

Wraps fastp with adapter auto-detection and handles paired-end and
single-end libraries automatically. Writes trimmed reads to
`data/reads_trimmed/`.

---

### Step 5 — Assemble reads into contigs

```bash
paleoamp assemble run \
    --reads-dir data/reads/ \
    --qc-report results/qc/qc_report.json \
    --output-dir results/assembly/ \
    --threads 8
```

Samples that failed aDNA QC are automatically skipped.

---

### Step 6 — Predict open reading frames

```bash
paleoamp predict orfs \
    --assembly-dir results/assembly/ \
    --output-dir results/orfs/
```

Runs Prodigal in metagenome mode. Outputs `.faa` (protein FASTA) and `.gff`
(gene coordinates) per sample.

---

### Step 7 — Amphipathicity pre-filter

Retains ORFs with an AMP-like amphipathic profile before the ML step.

```bash
paleoamp predict hydrophobic \
    --orfs-dir results/orfs/ \
    --output-dir results/orfs_hydrophobic/
```

Uses the **Eisenberg hydrophobic moment (µH)** — the vector sum of
per-residue hydrophobicities around a helical wheel at 100°/residue —
rather than mean hydrophobicity. µH captures the spatial segregation of
hydrophobic and hydrophilic faces that defines amphipathic AMPs, and is
insensitive to the cancellation that occurs when averaging both faces.

Default thresholds: µH ≥ 0.10, hydrophobic residue fraction ≥ 30%.

Note: this filter is tuned for α-helical AMPs. β-sheet AMPs (defensins,
protegrins) may not pass the µH threshold and will be missed.

---

### Step 8 — Score ORFs with the AMP classifier

```bash
paleoamp ml score \
    --orfs-dir results/orfs_hydrophobic/ \
    --checkpoint results/ml/checkpoints/best_model.pt \
    --output-dir results/ml/scores/
```

Generates frozen ESM-2 embeddings (480-dim) per sequence and scores them
with a trained MLP (480 → 256 → 64 → 1). Use `--all-samples` to score all
samples in one pass.

---

### Step 9 — Screen for novelty

```bash
paleoamp ml novelty \
    --amp-db data/amp_databases/merged_amp_db.fasta \
    --output-dir results/ml/novelty/
```

Compares high-scoring candidates against the merged AMP database using
MMseqs2 at `-s 7.5` sensitivity (full Smith-Waterman, no k-mer prefilter).
The prefilter is disabled because short peptides (10–50 aa) have too few
6-mers to shortlist reliably, producing false "novel" calls at default
sensitivity.

Sequences with no hit at ≥ 40% identity and ≥ 80% query coverage are
reported as candidate novel ancient AMPs.

---

### Step 10 — Validate candidates

```bash
paleoamp ml validate \
    --gff-dir results/orfs/ \
    --dataset data/ml/dataset_split.csv \
    --output-dir results/ml/validated/
```

Three independent validation layers:

1. **Physicochemical filters** — net charge ≥ 0, µH ≥ 0.10, length 10–100 aa,
   Shannon entropy ≥ 1.5, no simple repeats, all standard amino acids.
2. **Prodigal GFF partial-ORF scoring** — complete ORF (both ends) scores 1.0;
   one-end truncated (contig edge) scores 0.5; fragment (no start or stop) scores
   0.0. In aDNA assemblies, ~98% of ORFs are edge-truncated, so a binary
   complete/partial flag would eliminate nearly all candidates.
3. **AAC logistic regression** — a second classifier trained on amino acid
   composition, fully independent of the ESM-2 model (threshold 0.35).
4. **Boman index** — threshold calibrated from the 10th percentile of training
   positives (~−0.35). The hardcoded value of 0.5 used in earlier versions was
   calibrated against transmembrane helices and incorrectly excluded almost all
   amphipathic AMP candidates.

**Scoring:** each of the five criteria contributes 0.0–1.0 (ORF completeness
is the only non-binary one). **PASS** requires ≥ 4.0, **WARN** ≥ 2.0, **FAIL** < 2.0.
The output TSV includes a `failed_criteria` column showing exactly which
criteria each WARN/FAIL candidate missed.

---

## Classifier Training

Run these commands once to train the ESM-2 + MLP classifier.

```bash
# 1. Download training sequences
#    Positives → data/amp_sources/  (shared with novelty reference DB)
#    Negatives → data/ml/raw/       (non-AMP, transmembrane, signal, shuffled)
paleoamp ml collect

# 2. Cluster at 40% identity and split train/test without data leakage
paleoamp ml cluster --dataset data/ml/dataset.csv

# 3. Generate frozen ESM-2 embeddings (GPU recommended)
paleoamp ml embed --dataset data/ml/dataset_split.csv

# 4. Train the MLP (< 50 epochs with early stopping; val AUPR typically ~0.986)
paleoamp ml train --embeddings-dir data/ml/embeddings/ --output-dir results/ml/
```

**Positive sources** (label = 1):
- APD2024a — natural AMPs
- DRAMP 3.0 general and specific — experimentally validated
- UniProt KW-0929 — reviewed AMPs (all, bacterial, and archaeal subsets)
- BACTIBASE — experimentally validated bacteriocins

**Negative sources** (label = 0):
- UniProt non-AMP reviewed peptides
- UniProt short transmembrane proteins (hard negatives — share hydrophobic
  character with AMPs but lack cationic amphipathicity)
- UniProt short signal peptides (hard negatives — amphipathic but not antimicrobial)
- Shuffled positives (same AA composition as real AMPs, random order — forces
  the model to learn positional/structural patterns)

---

## Control Validation

To verify the pipeline is correctly calibrated, run known AMPs and shuffled
controls through the full pipeline:

```bash
# Controls FASTA is provided in data/controls/controls.faa
paleoamp predict hydrophobic -i results/orfs/controls -o results/orfs_hydrophobic
paleoamp ml score --all-samples --output-dir results/ml/scores
paleoamp ml novelty --scores-tsv results/ml/scores/controls_amp_scores.tsv
```

**Expected behaviour:**

| Control | Expected | Observed (val AUPR 0.986 model) |
|---|---|---|
| Magainin-2 | score ≥ 0.5, flagged known | 0.978, known ✓ |
| LL-37 | score ≥ 0.5, flagged known | 1.000, known ✓ |
| Melittin | score ≥ 0.5, flagged known | 1.000, known ✓ |
| Temporin-L | score ≥ 0.5, flagged known | 0.939, known ✓ |
| Cecropin-A | score ≥ 0.5, flagged known | 0.999, known ✓ |
| HNP-1 defensin | fails µH filter (β-sheet) | filtered ✓ |
| Shuffled magainin-2 | score < 0.5 | 0.004 ✓ |
| Shuffled temporin-L | score < 0.5 | 0.225 ✓ |
| Shuffled cecropin-A | score < 0.5 (known limitation) | **0.944 ✗** |

The cecropin-A shuffled false positive highlights a known limitation: sequences
with very AMP-like amino acid composition (high K/L/W) score high even when
disordered. Structural prediction (ESMFold) would distinguish these cases; the
AAC classifier partially compensates (cecropin-A shuffled fails AAC at 0.28).

---

## Configuration

Default thresholds live in `config/defaults.yaml`.

| Section | Parameter | Default | Description |
|---|---|---|---|
| `quality` | `min_mean_phred` | 20 | Minimum mean base quality |
| `quality` | `min_read_length` | 30 bp | Minimum read length |
| `adna_damage` | `min_ct_rate_5prime` | 0.10 | 5′ C→T damage threshold |
| `adna_damage` | `min_ga_rate_3prime` | 0.10 | 3′ G→A damage threshold |
| `assembly` | `min_contig_len` | 100 bp | Minimum assembled contig length |

---

## Output Files

| Path | Contents |
|---|---|
| `results/qc/qc_report.json` | Per-sample aDNA QC verdicts |
| `results/assembly/<id>/final.contigs.fa` | Assembled contigs |
| `results/orfs/<id>/<id>.faa` | Predicted proteins |
| `results/orfs/<id>/<id>.gff` | Gene coordinates (partial-ORF flags) |
| `results/orfs_hydrophobic/<id>/<id>.faa` | Amphipathicity-filtered proteins |
| `results/ml/scores/<id>_amp_scores.tsv` | AMP probability per ORF |
| `results/ml/novelty/<id>_novel.tsv` | Candidate novel AMPs (no DB hit) |
| `results/ml/novelty/<id>_known_hits.tsv` | ORFs matched to known AMPs |
| `results/ml/validated/<id>_novel_validated.tsv` | Validated candidates with PASS/WARN/FAIL verdict and `failed_criteria` annotation |
| `data/amp_sources/` | Downloaded AMP source FASTAs (shared by DB and ML) |
| `data/amp_databases/merged_amp_db.fasta` | Merged novelty-screen reference |
| `data/ml/dataset_split.csv` | Labelled train/test dataset |
| `results/ml/checkpoints/best_model.pt` | Trained MLP weights |
