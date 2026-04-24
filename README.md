# sge-utils

A collection of tools for Saturation Genome Editing (SGE) workflows, covering oligo library design and fitness score visualization.

---

## Repository Structure

```
sge-utils/
├── notebook_utils/           # Jupyter notebooks for SGE library design
├── SimpleSGEViz/             # CLI pipeline for SGE fitness score visualization
└── useful_scripts/           # Standalone utilities for SGE data processing and visualization
```

---

## Setup

A conda environment file is provided with all dependencies:

```bash
conda env create -f conda_sge-utils.yml
conda activate sge-utils
```

For the visualization pipeline, also install the `sgeviz` package:

```bash
cd SimpleSGEViz
pip install -e .
```

---

## notebook_utils

Interactive Jupyter notebooks for generating SGE oligo libraries and clonal homology arms for golden gate cloning.

### SGE_ClonalHA_Generator

Generates clonal homology arm (HA) sequences and Golden Gate AMP primers from WT homology arm sequences, for cloning SGE libraries via Golden Gate assembly

**Inputs:** A multi-sequence FASTA with 3 entries per target — forward AMP primer, reverse AMP primer, and full homology arm sequences, named `{GENE}_{TARGET}_AMP_F`, `{GENE}_{TARGET}_AMP_R`, and `{GENE}_{TARGET}_HA`.

**What it does:**
1. Pre-checks all sequences for restriction enzyme (REase) sites
2. Trims the full HA to AMP primer binding sites to produce 5' and 3' HA fragments
3. Builds clonal HA: `5' HA + stuffer sequence (with PaqCI sites) + 3' HA`
4. Builds Golden Gate AMP primers with standard GG stuffer overhangs
5. Validates outputs and saves clonal HAs and GG AMP primers as separate FASTA files

> The notebook uses PaqCI (CACCTGC) by default. Adjust the REase configuration if using a different enzyme.

---

### SGEoligos_SNVandFullMut_libs

Generates oligonucleotide libraries for SGE experiments, producing two complementary libraries:
- **SNV library** — all single-nucleotide variants and 3-bp deletions
- **Full saturation mutagenesis library** — all possible amino acid substitutions at every codon (optional output)

**Inputs:** An Excel (`.xlsx`) file with columns:
- `Oligo Name` — following the `{GENE}_{REGION}_*` convention
- `Oligo Sequence` — DNA sequence (lowercase bases = fixed edits such as PAM edits or HDR markers, which are preserved but excluded from mutagenesis)
- `Type` (optional) — `AMP Forward`, `AMP Reverse`, or `SGEoligo`
- `Library type` — `SNVs+3bpdel` or `SNVs+3bpdel+Saturation`

**What it does:**
1. Generates all SNVs and 3-bp deletions for each oligo
2. For saturation regions, fetches the canonical CDS from the Ensembl REST API, auto-detects reading frame and AA offset, and generates all missense variants, synonymous variants, stop codons, and 3-nt deletions
3. Handles intronic flanking sequences automatically
4. Outputs oligo libraries as text files with name and sequence

> Saturation mutagenesis requires network access to the Ensembl REST API.

---

## SimpleSGEViz

A command-line pipeline for generating standard visualization figures from SGE fitness score data. See [SimpleSGEViz/README.md](SimpleSGEViz/README.md) for full documentation.

### Quick start

```bash
sgeviz <input_dir> <output_dir> [--format html|png|svg] [--excel]
```

**Example:**

```bash
sgeviz ./data/BRCA1/ ./output/BRCA1/ --format png --excel
```

### What it generates

The pipeline auto-detects genes in the input directory and produces figures for each one. Core figures are always generated; additional figures are produced automatically when optional input files are present.

| Figure | Condition |
|---|---|
| Score distribution histogram + strip plot | Always |
| Replicate correlation heatmap | Always |
| Scores across gene (per-exon scatter) | Always |
| Exon or library cartoon | Always (requires Ensembl access) |
| Amino acid substitution heatmap | If `amino_acid_change` column present |
| Pathogenicity predictor scatter plots | If AlphaMissense / REVEL / CADD / MutPred2 scores present |
| ClinVar strip plot + ROC curve | If ClinVar file detected |
| Allele frequency vs score heatmap | If gnomAD or Regeneron file detected |
| Edit rate bar chart | If edit rates file detected |

### Required inputs


| File pattern | Description |
|---|---|
| `*{gene}allscores.tsv` | Combined SNV and 3bp deletion fitness scores |
| `*{gene}modelparams.tsv` | SGE model thresholds |
| `*{gene}snvcounts.tsv` | Per-replicate SNV read counts |
| `*{gene}delcounts.tsv` | Per-replicate 3bp deletion read counts |

Optional files (ClinVar, gnomAD, domain annotations, VEP output, edit rates, targets) are auto-detected by filename. See [SimpleSGEViz/README.md](SimpleSGEViz/README.md) for full input/output specifications.

---

## useful_scripts

Standalone scripts for common SGE data processing and structure visualization tasks.

### SGEColor_ChimeraX_MissenseOnly.py

Colors a protein ribbon structure in [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) by per-residue SGE scores using only missense variants. Scores are aggregated per amino acid position and mapped onto a white → red color scale clamped to [−0.2, 0].

> **Requirement:** [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) must be installed separately — it is not part of the conda environment. Required Python packages (`pandas`, `openpyxl`, `matplotlib`) are installed automatically into ChimeraX's Python environment on first run.

**Usage:**

In the ChimeraX command line:
```
runscript /path/to/SGEColor_ChimeraX_MissenseOnly.py
```

The script guides the user through a series of dialogs at runtime:

| Step | Dialog | Notes |
|---|---|---|
| 1 | PDB ID | Skipped if a structure is already loaded |
| 2 | SGE score file | File picker; accepts `.xlsx`, `.tsv`, or `.csv` |
| 3 | Chain selection | Dropdown populated from chains in the loaded structure |
| 4 | RNA score filter | Optional — exclude variants below a given `RNA_score` threshold |
| 5 | Add another? | Repeat steps 2–4 to color multiple chains in one run |

**Score file format:**

Required columns: `variant_qc_flag`, `consequence`, `amino_acid_change`, `score`

Optional column: `RNA_score` (used by the RNA score filter)

For `.xlsx` files, scores must be on a sheet named `scores`.

**Configuration** (edit at top of script):

| Parameter | Options | Default | Description |
|---|---|---|---|
| `analysis_type` | `'med'` \| `'mean'` \| `'min'` | `'med'` | Score aggregation method per residue |
| `show_legend` | `True` \| `False` | `False` | Show colorbar legend window |
| `save_legend` | `True` \| `False` | `False` | Prompt to save legend as PNG |
| `dna_style` | `'stubs'` \| `'slab'` \| `'fill'` \| `'atoms'` | `'stubs'` | ssDNA display style |

**Surface coloring:** The script targets atoms, cartoon, and surface (`target abcs`). For surface coloring to apply, show the surface before running the script (`surface` in the ChimeraX command line, or **Molecule Display → Surfaces → Show**).

---

### GetEditRates.sh

Computes per-target edit rates from readstats TSV files and writes them to a single output TSV.

**Usage:**

```bash
./GetEditRates.sh -i <input_directory> -o <output_directory>
```

**Arguments:**

| Flag | Description |
|---|---|
| `-i` | Path to the input directory (expects `*_X*_R*_D05*readstats.tsv` files inside) |
| `-o` | Path to the output directory where the edit rate file will be written |

The gene name and sample date are parsed from the input path. The output file is named `{gene}.editrates.{date}.tsv`.

**Output format:**

A tab-separated file with columns:
- `target_rep` — target/replicate identifier (column 1 from the readstats file)
- `edit_rate` — computed as `(col10 + col11) / col3` from each readstats file
