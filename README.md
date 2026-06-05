# sge-utils

A collection of tools for Saturation Genome Editing (SGE) workflows, covering oligo library design and fitness score visualization.

---

## Repository Structure

```
sge-utils/
├── examples/
│   ├── example_inputs/       # Example input files for the pipeline and notebooks
│   └── example_outputs/      # Expected outputs from the example inputs
├── notebook_utils/
│   ├── standalone_sgeviz/    # Standalone notebooks for individual figure generation
│   └── ...                   # SGE library design notebooks
├── SimpleSGEViz/             # CLI pipeline for SGE fitness score visualization
└── useful_scripts/           # Standalone utilities for SGE data processing and visualization
```

---

## Setup

**1. Clone the repository**

```bash
git clone https://github.com/ivanw314/sge-utils.git
cd sge-utils
```

**2. Install conda**

If conda is not already installed, download and install [Miniconda](https://docs.anaconda.com/miniconda/) (recommended) or [Anaconda](https://www.anaconda.com/download).

**3. Create and activate the environment**

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

## Examples

Example input and output files for BARD1 are provided in `examples/` to help you verify your setup and explore expected outputs.

### sgeviz pipeline

Inputs in `examples/example_inputs/sgeviz/`:

| File | Description |
|---|---|
| `BARD1.allscores.v1.2.1.tsv` | Combined SNV and 3bp deletion fitness scores |
| `BARD1.modelparams.v1.2.1.tsv` | Fitness score thresholds |
| `BARD1.snvcounts.tsv` | Per-replicate SNV read counts |
| `BARD1.delcounts.tsv` | Per-replicate 3bp deletion read counts |
| `BARD1.editrates.tsv` | Per-target edit rates |
| `BARD1.targets.tsv` | Coordinates for library amplicons |
| `BARD1_domains.xlsx` | Protein domain annotations |
| `20240905_BARD1_gnomADv4.1.0_SNVs.xlsx` | gnomAD allele frequencies |
| `20250912_BARD1_ClinVarSNVs_1StarPlus.txt` | ClinVar variant classifications |
| `20251002_BARD1snvs_VEP.xlsx` | VEP output with computational predictor scores |

To run the pipeline on the example data (working directory should be the repo root, `sge-utils/`):

```bash
sgeviz examples/example_inputs/sgeviz/ examples/example_outputs/sgeviz/ --format png --excel
```

Expected outputs are in `examples/example_outputs/sgeviz/` — `BARD1_data.xlsx` plus 13 figure PNGs in `fig_outputs/`.

### SGE_ClonalHA_Generator

- Input: `examples/example_inputs/clonalHA_generation/20260316_HA_clonalDNA_subset_HA1.fa`
- Outputs: `examples/example_outputs/clonalHA_generation/` — clonal HA sequences and Golden Gate AMP primers as FASTA files

### SGEoligos_SNVandFullMut_libs

- Input: `examples/example_inputs/library_generation/20260317_SGE_All_Oligos_merged_recoded.xlsx`
- Output: `examples/example_outputs/library_generation/` — oligo library as a text file

---

## notebook_utils

### standalone_sgeviz

Standalone Jupyter notebooks for generating individual `sgeviz` figures without running the full pipeline. Each notebook loads data directly from the multi-sheet Excel workbook produced by `sgeviz --excel` (`*_data.xlsx`) and saves a PNG by default (change the extension to `.html` for interactive output or `.svg` as needed).

| Notebook | Figure | Required sheets |
|---|---|---|
| `histogram_stripplot.ipynb` | Score distribution histogram + strip plot | `scores`, `thresholds` |
| `correlation_heatmap.ipynb` | Replicate Pearson r heatmap | `counts` |
| `scores_across_gene.ipynb` | Per-exon score scatter | `scores`, `thresholds` |
| `aa_heatmap.ipynb` | Amino acid substitution heatmap | `scores`, `thresholds` |
| `predictor_scatter.ipynb` | SGE vs. computational predictor scatter | `scores`, `thresholds` |
| `clinvar.ipynb` | ClinVar strip plot + ROC curve | `scores`, `thresholds` |
| `rna_score.ipynb` | RNA score scatter + stem plot | `scores`, `thresholds` |
| `edit_rate_barplot.ipynb` | Library edit rate barplot | `edit_rates` |
| `maf_vs_score.ipynb` | Allele frequency vs. fitness score heatmap | `scores` (requires `gnomad_af`/`regeneron_af` columns) |
| `gene_cartoon.ipynb` | Gene exon structure or library amplicon cartoon | Ensembl API (internet required) |
| `single_track_cartoon.ipynb` | Single-track library cartoon (matplotlib) | Ensembl API + targets TSV |

Each notebook has a `# --- Plot customization (optional) ---` block in its configuration cell exposing width, height, and relevant axis/color domain parameters (e.g. `score_domain`, `rna_domain`, `panel_width`). These map directly to new keyword arguments on the underlying figure functions and default to the same values used by the pipeline.

**Requirements:** `sge-utils` conda environment with `sgeviz` installed. PNG/SVG output additionally requires `vl-convert-python` (`pip install vl-convert-python`).

#### Usage

**VS Code**

1. Open the notebook file (`.ipynb`) in VS Code.
2. When prompted to select a kernel, choose the `sge-utils` conda environment. If it doesn't appear, install the [Jupyter extension](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) and ensure the environment is registered with `python -m ipykernel install --user --name sge-utils`.
3. Edit the `# --- Configuration ---` cell: set `excel_path` (or `gene`/`targets_path` for cartoon notebooks) to point to your files.
4. Run all cells with **Run All** (⇧⌘↩ on Mac, Ctrl+Alt+Enter on Windows/Linux) or step through cells individually with ⇧↩.

**JupyterLab**

1. Activate the environment and launch JupyterLab:
   ```bash
   conda activate sge-utils
   jupyter lab
   ```
2. Navigate to `notebook_utils/standalone_sgeviz/` in the file browser and open the desired notebook.
3. Edit the `# --- Configuration ---` cell with your file paths, then run all cells via **Run → Run All Cells**.

---

### Library design notebooks

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

### Figure customization

All figure functions in `sgeviz.figures` now accept explicit `width`, `height`, and domain parameters so that individual plots can be resized or re-ranged without touching the pipeline. The pipeline itself continues to use the original defaults. The standalone notebooks in `notebook_utils/standalone_sgeviz/` expose these parameters in a dedicated configuration block.

| Parameter | Applies to | Default |
|---|---|---|
| `width`, `height` | All figures | varies per figure |
| `score_domain` | `aa_heatmap`, `rna_score.make_scatter` | `(-0.2, 0)` / `(-0.6, 0.3)` |
| `rna_domain` | `rna_score.make_scatter`, `rna_score.make_stem_plot` | `(-8, 3)` |
| `height_per_row` | `aa_heatmap` | `25` |
| `panel_width`, `panel_height` | `predictor_scatter` | `350`, `250` |
| `bar_width` | `edit_rate_barplot` | `35` |
| `width`, `height` | `maf_score.make_plot` | `300`, `250` |
| `width` | `gene_cartoon.make_exon_cartoon`, `make_library_cartoon` | `800` |
| `fig_width` | `single_track_cartoon.make_plot` | `10` |

---

## useful_scripts

Standalone scripts for common SGE data processing and structure visualization tasks.

### SGEColor_ChimeraX_MissenseOnly.py

Colors a protein ribbon structure in [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) by per-residue SGE scores using only missense variants. Scores are aggregated per amino acid position and mapped onto a configurable white ↔ red color scale. The script is fully dialog-driven — no command-line arguments needed.

> **Requirement:** [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) must be installed separately — it is not part of the conda environment. Required Python packages (`pandas`, `openpyxl`, `matplotlib`) are installed automatically into ChimeraX's Python environment on first run.

**Usage:**

In the ChimeraX command line:
```
runscript /path/to/SGEColor_ChimeraX_MissenseOnly.py
```

**Dialog sequence:**

| Step | Dialog | Notes |
|---|---|---|
| 1 | PDB ID | Skipped if a structure is already open |
| 2 | Aggregation method | Median / mean / min per-residue aggregation |
| 3 | Score settings | Default or custom: color range, color direction (score column deferred to file-load in custom mode) |
| 4 | Colorbar legend | Don't show / show / show and save |
| 5 | SGE score file | File picker; accepts `.xlsx`, `.tsv`, `.csv` |
| — | Sheet selection | Shown for multi-sheet `.xlsx` files; `scores` pre-selected if present |
| — | Score column | Shown in custom mode after file/sheet load; dropdown of all columns in the file |
| 6 | Chain selection | Dropdown of chains in the loaded structure |
| 7 | RNA score filter | Optional — cancel to skip |
| 8 | Add another? | Repeat steps 5–7 for additional chains |
| — | Column remapping | Appears per file if expected column names are not found |
| — | Residue range check | Confirm matched range before coloring proceeds |
| — | Offset detection | Shown if sequence identity < 90% and a better offset exists |
| — | Validation gate | Skip or force-proceed if identity < 80% after offset check |

**Score file format:**

Accepted file types: `.xlsx`, `.tsv`, `.csv`. For `.xlsx` files with multiple sheets a picker dialog appears; the sheet named `scores` is pre-selected if present.

| Column | Required | Description |
|---|---|---|
| `consequence` | Yes | Rows containing `'missense_variant'` are kept; other rows ignored |
| `amino_acid_change` | Yes | Amino acid substitution — accepts `A123G`, `p.Met123Val`, `NP_xxx:p.Met123Lys`, and other standard HGVS forms |
| `score` | Yes (default name) | Numeric score used for coloring; name is configurable |
| `variant_qc_flag` | No | Rows where this equals `'WARN'` are excluded if the column exists |
| `RNA_score` / `RNAscore` | No | Used by the optional RNA score filter |

If any required column name differs in your dataset, a picker dialog appears automatically showing all available columns so you can select the correct one.

**Validation:**

For each chain, the script checks coverage (how many scored positions are found in the structure) and sequence identity (whether the reference amino acid in the score file matches the PDB residue at each position). A residue range confirmation dialog is shown before coloring. If identity is imperfect, an automatic scan for a residue number offset is run (common when the PDB starts at a different residue than the canonical sequence, e.g. after signal peptide removal). Non-standard residues in the structure (e.g. MSE, SEP, TPO, HYP) are mapped to their canonical parent amino acid before comparison so they do not inflate the mismatch count. Positions with amino acid mismatches are always left gray. Coloring is blocked by default if identity is below 80% and no offset corrects it. A validation summary is printed to the log after coloring finishes.

**Configuration** (edit at top of script; all settings are also available via the GUI):

| Parameter | Options | Default | Description |
|---|---|---|---|
| `analysis_type` | `'med'` \| `'mean'` \| `'min'` | `'med'` | Score aggregation method per residue |
| `score_column` | str | `'score'` | Column name in the score file to use |
| `clamp_min` | float | `-0.2` | Lower bound of the color normalization range |
| `clamp_max` | float | `0.0` | Upper bound of the color normalization range |
| `high_is_red` | `True` \| `False` | `False` | If `False`, white = high score, red = low; if `True`, reversed |
| `show_legend` | `True` \| `False` | `False` | Show colorbar legend window |
| `save_legend` | `True` \| `False` | `False` | Prompt to save legend as PNG |
| `dna_style` | `'stubs'` \| `'slab'` \| `'fill'` \| `'atoms'` | `'stubs'` | ssDNA display style |

**Surface coloring:** The script targets atoms, cartoon, and surface (`target abcs`). For surface coloring to apply, show the surface **before** running the script (`surface` in the ChimeraX command line, or **Molecule Display → Surfaces → Show**).

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
