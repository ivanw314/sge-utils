"""
SGEColor_ChimeraX_MissenseOnly.py
==================================
Colors a protein ribbon structure in ChimeraX by per-residue SGE (Saturating Genome Editing)
scores, using only missense variants. Scores are aggregated per amino acid position and mapped
onto a color scale (white ↔ red) whose range and direction are configurable.

USAGE
-----
In the ChimeraX command line:
    runscript /path/to/SGEColor_ChimeraX_MissenseOnly.py

INPUT FILE FORMAT
------------------
Supported file types: .xlsx, .tsv, .csv
  For .xlsx: if the file has multiple sheets, a picker dialog is shown;
             the sheet named 'scores' is pre-selected if present.

Required columns:
  consequence        str    Only rows containing 'missense_variant' are kept.
  amino_acid_change  str    Amino acid substitution in any of these formats:
                              A123G             (one-letter shorthand)
                              p.A123G           (HGVS one-letter)
                              p.Met123Val       (HGVS three-letter)
                              p.Met123V         (HGVS mixed)
                              NP_009225.1:p.Met123Lys  (accession-prefixed HGVS)
                            The ref AA is used to verify the sequence matches
                            the selected PDB chain. Non-standard residues in
                            the structure (e.g. MSE, SEP, TPO) are mapped to
                            their canonical equivalents before comparison.
  <score_column>     float  Numeric score used for coloring. Column name is set by
                            the score_column config variable (default: 'score').

Optional quality-control columns:
  variant_qc_flag    str    If present, rows where this equals 'WARN' are excluded.

Optional columns:
  RNA_score / RNAscore  float  If an RNA score threshold is specified at runtime,
                               variants below the threshold are excluded.
                               Rows with a missing RNA score are kept regardless.

INTERACTIVE DIALOGS (shown at runtime)
---------------------------------------
Shown once per run:
  1. PDB ID           — prompted only if no structure is currently loaded.
                        If a model is already open, that structure is used as-is.

Repeated per gene/chain (steps 2–7 repeat for each additional chain):
  2. Aggregation      — choose median / mean / min per-residue aggregation.
  3. Score settings   — use defaults or override color range and direction.
                        If Custom is chosen, score column is deferred to step 5
                        so it can be picked from the actual file's columns.
  4. Colorbar legend  — don't show / show / show and save.
  5. SGE score file   — file picker; accepts .xlsx, .tsv, .csv.
     Sheet selection  — shown immediately after for multi-sheet .xlsx files.
     Score column     — shown immediately after (custom mode only): dropdown of
                        all columns in the loaded file.
  6. Chain selection  — dropdown of chains in the loaded structure.
  7. RNA filter       — optional threshold; cancel to skip.
  8. Add another?     — if yes, loop restarts from step 2.

Shown only when columns are missing:
  - Column Not Found  — picker over all available columns, shown separately
                        for consequence, amino_acid_change, and score_column.
  - Missense label    — picker over unique consequence values if no
                        'missense_variant' rows are found after column mapping.
  - RNA column        — same picker if RNA_score is not found and a threshold
                        was requested.

Shown per chain during validation:
  9. Residue range    — confirm matched residue range before coloring proceeds.
  10. Offset detected — apply a detected residue number offset (e.g. signal
                        peptide removed), if one is found. Conditional.
  11. Validation gate — skip or proceed if sequence identity is below 80%.
                        Conditional; default is skip.

SURFACE COLORING
-----------------
The script automatically generates a molecular surface for each colored chain
before applying scores. To hide the surface and show only the ribbon diagram,
run in the ChimeraX command line or use Molecule Display → Surfaces → Hide:
    hide /{chain_id} surface

VALIDATION
----------
For each chain, the script checks:
  - Coverage: how many scored positions are found in the chain, and the
    matched residue range (shown as a confirmation dialog).
  - Sequence identity: whether the reference amino acid in the score file
    matches the structure at each position. Mismatched positions are left
    gray even if the user proceeds. If identity is below 90%, an automatic
    scan for a residue number offset is run.
  - Coloring is blocked (with an override option) if identity is below 80%
    and no offset corrects it.
A validation summary is printed to the log after coloring completes.

CONFIGURATION (edit at top of script — all can also be set via the GUI)
---------------------------------------
  analysis_type  'med' | 'mean' | 'min'   Aggregation method per residue (default: 'med')
  score_column   str                       Column name in the score file to use (default: 'score')
  clamp_min      float                     Lower bound of the color normalization range (default: -0.2)
  clamp_max      float                     Upper bound of the color normalization range (default: 0)
  high_is_red    True | False              If True, high scores → red, low → white.
                                           If False (default), high scores → white, low → red.
  show_legend    True | False              Show colorbar legend window (default: False)
  save_legend    True | False              Prompt to save legend as PNG (default: False)
  dna_style      'stubs'|'slab'|'fill'    ssDNA display style (default: 'stubs')
                 |'atoms'
"""

from chimerax.core.commands import run
import subprocess
import sys
import os
import re
from Qt.QtWidgets import QInputDialog, QFileDialog

#Install required packages into ChimeraX's Python environment if not already present
for _pkg in ['pandas', 'openpyxl', 'matplotlib']:
    try:
        __import__(_pkg)
    except ImportError:
        print(f'Installing {_pkg} into ChimeraX Python environment...')
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', _pkg])

import pandas as pd
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

analysis_type      = 'med'   # 'med', 'mean', or 'min' — aggregation method per residue
score_column       = 'score' # column in the score file to use for coloring
custom_score_mode  = False   # True when the user chose Custom settings; triggers score column picker at file-load time
clamp_min     = -0.2    # lower bound for color normalization range
clamp_max     = 0.0     # upper bound for color normalization range
high_is_red   = False   # False → high score = white, low = red (default); True → high = red, low = white
show_legend = False     # whether to show the legend figure
save_legend = False     # whether to save the legend figure
dna_style = 'stubs'     # ssDNA display style: 'stubs', 'slab', 'fill', or 'atoms'


def resolve_column(df, expected_name, description, parent, force_pick=False):
    """Return expected_name if present (and force_pick is False); otherwise show a picker so the user can select the right column."""
    if not force_pick and expected_name in df.columns:
        return expected_name
    if force_pick:
        title = f'Select Score Column'
        msg   = f'Select the column to use as the score ({description}):'
        default_idx = df.columns.tolist().index(expected_name) if expected_name in df.columns else 0
    else:
        title = f'Column Not Found: "{expected_name}"'
        msg   = (f'Expected column "{expected_name}" ({description})\n'
                 f'was not found in this file.\n\n'
                 f'Select the column that corresponds to it:')
        default_idx = 0
    chosen, ok = QInputDialog.getItem(parent, title, msg, df.columns.tolist(), default_idx, False)
    if not ok:
        raise ValueError(f'Required column "{expected_name}" not found and not remapped — aborting.')
    return chosen


def read_scores(file, parent, rna_score_threshold=None): #Reads score file
    global score_column, custom_score_mode
    ext = os.path.splitext(file)[1].lower()
    if ext == '.xlsx':
        xl = pd.ExcelFile(file)
        sheets = xl.sheet_names
        if len(sheets) == 1:
            sheet = sheets[0]
        else:
            default_idx = sheets.index('scores') if 'scores' in sheets else 0
            sheet, ok = QInputDialog.getItem(
                parent, 'Select Sheet',
                f'"{os.path.basename(file)}" has multiple sheets.\nSelect the one containing SGE scores:',
                sheets, default_idx, False)
            if not ok:
                raise ValueError('Sheet selection cancelled — aborting.')
        df = pd.read_excel(xl, sheet_name=sheet)
    elif ext == '.tsv':
        df = pd.read_csv(file, sep='\t')
    elif ext == '.csv':
        df = pd.read_csv(file)
    else:
        raise ValueError(f'Unsupported file type: {ext}. Use .xlsx, .tsv, or .csv')

    if 'variant_qc_flag' in df.columns:
        df = df.loc[df['variant_qc_flag'] != 'WARN'] #Filters out variants with WARN flag

    consequence_col = resolve_column(df, 'consequence',       'identifies missense variants', parent)
    aa_change_col   = resolve_column(df, 'amino_acid_change',
                                     'amino acid substitution — accepts A123G, p.Met123Val, NP_xxx:p.Met123Lys, etc.', parent)
    score_col       = resolve_column(df, score_column,        'numeric score used for coloring', parent,
                                     force_pick=custom_score_mode)
    if score_col != score_column:
        score_column = score_col  # keep legend label in sync

    df = df.rename(columns={consequence_col: 'Consequence', aa_change_col: 'AAsub', score_col: 'snv_score'})

    missense_filter = 'missense_variant'
    if not df['Consequence'].str.contains(missense_filter).any():
        unique_consequences = sorted(df['Consequence'].dropna().unique().tolist())
        chosen, ok = QInputDialog.getItem(
            parent,
            'Missense Label Not Found',
            f'No rows containing "{missense_filter}" were found in the consequence column.\n\n'
            f'Select how missense variants are labelled in your dataset:',
            unique_consequences, 0, False)
        if not ok:
            raise ValueError('Missense variant label not specified — aborting.')
        missense_filter = chosen

    df = df.loc[df['Consequence'].str.contains(missense_filter, regex=False)]

    if rna_score_threshold is not None: #Optionally filter out variants below the RNA score threshold; NaN rows are kept
        rna_col = next((c for c in df.columns if c in ('RNA_score', 'RNAscore')), None)
        if rna_col is None:
            rna_col = resolve_column(df, 'RNA_score', 'RNA score used for optional filtering', parent)
        df = df.loc[df[rna_col].isna() | (df[rna_col] >= rna_score_threshold)]

    parsed = df['AAsub'].map(parse_aa_change)
    n_unparseable = parsed.isna().sum()
    if n_unparseable:
        print(f'  Warning: {n_unparseable} row(s) with unparseable amino_acid_change values dropped.')
    df = df[parsed.notna()].copy()
    parsed = parsed[parsed.notna()]
    df['AApos']       = parsed.map(lambda x: x[1])
    df['ref_1letter'] = parsed.map(lambda x: x[0])

    return df

def group_scores(df): #Groups variant scores by AA position and creates calculates min and mean score

    df = df[['AApos', 'snv_score']] #Gets AA position and score column only

    grouped = df.groupby('AApos') #Groups by position

    min_scores = grouped['snv_score'].min().reset_index() #Summary dataframe with minimum scores
    mean_scores = grouped['snv_score'].mean().reset_index() #Summary dataframe with mean scores
    median_scores = grouped['snv_score'].median().reset_index() #Summary dataframe with median scores

    min_scores = min_scores.set_index('AApos')['snv_score'].to_dict() #Turns min scores dataframe into dictionary for coloring
    mean_scores = mean_scores.set_index('AApos')['snv_score'].to_dict() #Turns mean scores dataframe into dictionnary for coloring
    median_scores = median_scores.set_index('AApos')['snv_score'].to_dict() #Turns median scores dataframe into dictionnary for coloring

    return min_scores, mean_scores, median_scores


def normalize_values(values): #Normalizes all values between 0 and 1 for coloring
    # Clamp all values to [clamp_min, clamp_max]
    clamped_values = {k: min(max(v, clamp_min), clamp_max) for k, v in values.items()}
    clamped_values = {k: v for k, v in clamped_values.items() if not pd.isna(v)} #Filters out NA values

    # Avoid division by zero if all values are the same
    if clamp_max == clamp_min:
        return {k: 0 for k in clamped_values.keys()}

    # Perform normalization
    return {k: (v - clamp_min) / (clamp_max - clamp_min) for k, v in clamped_values.items()}


AA_ONE_TO_THREE = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
}

THREE_TO_ONE = {v: k for k, v in AA_ONE_TO_THREE.items()}

# Maps non-standard residue three-letter codes (as returned by ChimeraX r.name) to their
# canonical parent amino acid, so crystallographic modifications don't cause false mismatches.
NONSTANDARD_TO_CANONICAL = {
    'MSE': 'MET',  # selenomethionine — very common in X-ray structures
    'HYP': 'PRO',  # 4-hydroxyproline
    'SEP': 'SER',  # phosphoserine
    'TPO': 'THR',  # phosphothreonine
    'PTR': 'TYR',  # phosphotyrosine
    'CSO': 'CYS',  # S-hydroxycysteine
    'CSS': 'CYS',  # S-mercaptocysteine
    'CME': 'CYS',  # carboxymethyl cysteine
    'OCS': 'CYS',  # cysteic acid
    'MLY': 'LYS',  # N6-methyl-lysine
    'M3L': 'LYS',  # N6-trimethyl-lysine
    'ALY': 'LYS',  # N6-acetyl-lysine
    'KCX': 'LYS',  # N6-carboxylysine
    'HIC': 'HIS',  # 4-methyl-histidine
    'TYS': 'TYR',  # O-sulfo-tyrosine
    'NEP': 'HIS',  # N1-phosphohistidine
    'MHO': 'MET',  # S-oxymethionine
}

# Matches the following amino acid change formats (all produce ref_1letter + position):
#   A123G                     — one-letter shorthand
#   p.A123G                   — HGVS one-letter
#   p.Met123Val / p.Met123V   — HGVS three-letter ref (mixed alt ok)
#   NP_009225.1:p.Met123Lys   — accession-prefixed HGVS
_AA_SUB_RE = re.compile(
    r'^(?:[^:]+:)?'           # optional accession prefix (e.g. NP_009225.1:)
    r'(?:p\.)?'               # optional p. prefix
    r'([A-Za-z]{1,3})'        # ref AA (1 or 3 letters)
    r'(\d+)'                  # residue position
    r'([A-Za-z*=]{1,3})$'     # alt AA, stop (*), or synonymous (=)
)

def parse_aa_change(s):
    """Return (ref_1letter, pos_int) parsed from an amino acid change string, or None if unparseable."""
    if not isinstance(s, str):
        return None
    m = _AA_SUB_RE.match(s.strip())
    if not m:
        return None
    ref_raw, pos_str = m.group(1), m.group(2)
    if len(ref_raw) == 1:
        ref_1 = ref_raw.upper()
    else:
        ref_1 = THREE_TO_ONE.get(ref_raw.upper())
    if ref_1 not in AA_ONE_TO_THREE:
        return None
    return ref_1, int(pos_str)


colors = [(1.0, 1.0, 1.0), (1, 0, 0)]  # White -> Red
n_bins = 1000  # Number of bins for the colormap
cmap_name = 'gray_to_red'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins) #makes the color map

def create_colorbar_legend():
    plt.close('all')  # clear any figures left over from previous runs
    fig, ax = plt.subplots(figsize=(1, 0.5))
    fig.subplots_adjust(bottom=0.5)

    # Match the colormap direction used by get_color
    legend_cmap = custom_cmap if high_is_red else custom_cmap.reversed()

    norm = plt.Normalize(vmin=clamp_min, vmax=clamp_max)
    sm = cm.ScalarMappable(cmap=legend_cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, cax=ax, orientation='horizontal')
    cbar.set_ticks([clamp_min, clamp_max])
    cbar.set_label(score_column, labelpad=10)

    if show_legend:
        plt.show()
    return fig


def get_color(value): #Gets color for each residue from score
    # high_is_red=False (default): high value → white, low → red  (invert lookup)
    # high_is_red=True:            high value → red,   low → white (normal lookup)
    return custom_cmap(value if high_is_red else 1 - value)


def rgb_to_hex(r, g, b): #Converts float RGB (0-1) to hex string for ChimeraX
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


def find_best_offset(ref_aa_by_pos, chain_residue_map, max_offset=500):
    """Scan integer offsets [-max_offset, +max_offset] to find the shift that maximises
    AA identity between score-file positions and the chain. Returns (offset, n_matches, n_checkable)."""
    scored = [(pos, AA_ONE_TO_THREE[aa]) for pos, aa in ref_aa_by_pos.items()
              if aa in AA_ONE_TO_THREE]
    if not scored:
        return 0, 0, 0
    best_offset, best_matches = 0, 0
    for offset in range(-max_offset, max_offset + 1):
        matches = sum(1 for pos, exp in scored
                      if NONSTANDARD_TO_CANONICAL.get(
                          chain_residue_map.get(pos + offset, ''),
                          chain_residue_map.get(pos + offset, '')) == exp)
        if matches > best_matches:
            best_matches, best_offset = matches, offset
    return best_offset, best_matches, len(scored)


def get_score_config(session):
    """Ask aggregation method and legend preference (always); optionally override clamp range / color direction.
    When Custom is chosen, score column selection is deferred to file-load time so the user picks from actual columns."""
    global analysis_type, score_column, clamp_min, clamp_max, high_is_red, show_legend, save_legend, custom_score_mode
    custom_score_mode = False  # reset each run so Default doesn't inherit a previous Custom choice
    parent = session.ui.main_window

    # Aggregation method — always asked
    agg_choice, ok = QInputDialog.getItem(
        parent, 'Aggregation Method',
        'How should per-residue scores be aggregated across variants?',
        ['Median', 'Mean', 'Minimum'],
        ['med', 'mean', 'min'].index(analysis_type), False)
    if ok:
        analysis_type = {'Median': 'med', 'Mean': 'mean', 'Minimum': 'min'}[agg_choice]

    mode, ok = QInputDialog.getItem(
        parent, 'Score Settings',
        'Which score settings would you like to use?',
        [f'Default  (column: "{score_column}",  range: {clamp_min} to {clamp_max},  white = high score)',
         'Custom score column / range / color direction'],
        0, False)
    if ok and mode.startswith('Custom'):
        custom_score_mode = True
        # Score column is selected after the file is loaded so the user picks from actual column names.

        # Clamp range — min
        cmin, ok = QInputDialog.getDouble(
            parent, 'Color Range — Minimum',
            'Lower bound for color normalization\n(values below this are clamped to this color):',
            clamp_min, -1e6, 1e6, 4)
        if ok:
            clamp_min = cmin

        # Clamp range — max
        cmax, ok = QInputDialog.getDouble(
            parent, 'Color Range — Maximum',
            'Upper bound for color normalization\n(values above this are clamped to this color):',
            clamp_max, -1e6, 1e6, 4)
        if ok:
            clamp_max = cmax

        if clamp_min > clamp_max:
            clamp_min, clamp_max = clamp_max, clamp_min
            print(f'  Note: min > max — values swapped to [{clamp_min}, {clamp_max}].')

        # Color direction
        direction, ok = QInputDialog.getItem(
            parent, 'Color Direction',
            'Which end of the score range should be red?',
            ['White = high score  /  Red = low score',
             'Red = high score  /  White = low score'],
            0, False)
        if ok:
            high_is_red = direction.startswith('Red')

    # Legend — always asked
    default_legend_idx = 2 if save_legend else (1 if show_legend else 0)
    legend_choice, ok = QInputDialog.getItem(
        parent, 'Colorbar Legend',
        'Show a colorbar legend for the score range?',
        ["Don't show legend", 'Show legend', 'Show and save legend'],
        default_legend_idx, False)
    if ok:
        show_legend = legend_choice != "Don't show legend"
        save_legend = legend_choice == 'Show and save legend'


def get_gene_config(session, available_chains):
    """Collect file, chain, and RNA filter for a single chain. Returns (file_path, chain_id, rna_threshold) or None if cancelled."""
    parent = session.ui.main_window
    file_path, _ = QFileDialog.getOpenFileName(
        parent, 'Select SGE Score File', '',
        'Score Files (*.xlsx *.tsv *.csv);;Excel (*.xlsx);;TSV (*.tsv);;CSV (*.csv)')
    if not file_path:
        return None
    chain_id, ok = QInputDialog.getItem(
        parent, 'Select Chain',
        f'Chain for {os.path.basename(file_path)}:',
        available_chains, 0, False)
    if not ok:
        return None
    rna_threshold, ok = QInputDialog.getDouble(
        parent, 'RNA Score Filter (optional)',
        'Exclude variants where RNA_score is below (cancel to skip):',
        0.0, -10.0, 10.0, 3)
    return (file_path, chain_id, rna_threshold if ok else None)


def main():  # 'session' is injected as a global by ChimeraX at runtime via runscript
    parent = session.ui.main_window

    # Load structure first so the user has it in view while configuring score settings
    existing_models = session.models.list()
    if not existing_models:
        pdb_id, ok = QInputDialog.getText(parent, 'Load Structure', 'Enter PDB ID:')
        if not ok or not pdb_id.strip():
            raise ValueError('No PDB ID provided — script aborted.')
        pdb_id = pdb_id.strip().upper()
        print(f'Loading structure {pdb_id}...')
        run(session, f'open {pdb_id} from pdb') #Fetches PDB structure from RCSB

        run(session, 'show cartoons')  #Show cartoon for all chains in the model
        run(session, 'hide atoms')     #Hide atoms for all chains
        run(session, 'hide bonds')     #Hide bonds for all chains
        run(session, 'show ~ protein & ~ solvent & ~ nucleic atoms') #Show ligands, ATP, ions as sticks
        run(session, 'show ~ protein & ~ solvent & ~ nucleic bonds') #Show their bonds
        if dna_style == 'atoms':
            run(session, 'show nucleic atoms')  #Show ssDNA as atoms/sticks
            run(session, 'show nucleic bonds')
        else:
            run(session, f'nucleotides {dna_style}')  #Show ssDNA bases (stubs/slab/fill)
    else:
        print(f'Existing session detected ({len(existing_models)} model(s) open). Skipping structure load — applying colors only.')

    models = session.models.list()
    available_chains = sorted(set(c.chain_id for m in models if hasattr(m, 'chains') for c in m.chains))

    first_chain = True
    while True:
        get_score_config(session)

        legend = create_colorbar_legend()
        if save_legend:
            save_path, _ = QFileDialog.getSaveFileName(parent, 'Save Legend', '', 'PNG Files (*.png)')
            if save_path:
                legend.savefig(save_path, dpi=500)
        if not show_legend:
            plt.close(legend)

        config = get_gene_config(session, available_chains)
        if config is None:
            if first_chain:
                raise ValueError('No data files specified — script aborted.')
            break
        first_chain = False
        file_path, chain_id, rna_score_threshold = config

        label = os.path.basename(file_path)
        print(f'Reading SGE scores from {label}...')
        raw_scores = read_scores(file_path, parent, rna_score_threshold)
        ref_aa_by_pos = raw_scores.drop_duplicates('AApos').set_index('AApos')['ref_1letter'].to_dict()

        print('Grouping scores by residue...')
        min_scores, mean_scores, median_scores = group_scores(raw_scores) #Gets min/mean/median score dataframes

        if analysis_type == 'mean':
            residue_values = mean_scores
        elif analysis_type == 'min':
            residue_values = min_scores
        elif analysis_type == 'med':
            residue_values = median_scores
        else:
            raise ValueError("analysis_type must be 'mean', 'min', or 'med'")

        print('Normalizing values for coloring...')
        normalized_values = normalize_values(residue_values) #Scores normalized to between 0 and 1

        # Build residue map for the target chain from the loaded structure
        chain_residue_map = {}
        for m in models:
            if not hasattr(m, 'chains'):
                continue
            for c in m.chains:
                if c.chain_id == chain_id:
                    chain_residue_map = {r.number: r.name for r in c.residues if r is not None}
                    break
            if chain_residue_map:
                break

        # --- Validation ---
        # Informational messages are collected and printed after coloring so they
        # aren't buried under ChimeraX command output.
        validation_log = []
        serious_issue = False
        skip_chain = False
        n_scored  = len(normalized_values)
        common    = set(normalized_values) & set(chain_residue_map)
        missing   = sorted(set(normalized_values) - set(chain_residue_map))

        # Coverage — range dialog runs now (can abort); stats deferred to log
        coverage_pct = 100 * len(common) / n_scored if n_scored else 0
        validation_log.append(f'Coverage: {len(common)}/{n_scored} scored positions found in chain {chain_id} ({coverage_pct:.0f}%)')
        if common:
            range_ok, ok = QInputDialog.getItem(
                parent, 'Check Residue Range',
                f'Scored positions matched chain {chain_id} residues {min(common)}–{max(common)}\n'
                f'({len(common)} of {max(common) - min(common) + 1} residues in that range scored).\n\n'
                f'Does this range look correct?',
                ['Yes — continue', 'No — skip this chain'],
                0, False)
            if ok and range_ok.startswith('No'):
                print(f'Skipping coloring for chain {chain_id} — residue range rejected.')
                skip_chain = True
        if missing:
            preview = missing[:10]
            suffix  = f' ... and {len(missing) - 10} more' if len(missing) > 10 else ''
            validation_log.append(f'  Note: {len(missing)} scored position(s) not found in chain {chain_id}: {preview}{suffix}')

        if not skip_chain:
            # Amino-acid identity
            n_checked = n_mismatched = 0
            identity_pct = 0.0
            mismatch_positions = set()
            mismatch_details = []
            for pos in sorted(common):
                ref_1letter = ref_aa_by_pos.get(pos)
                expected = AA_ONE_TO_THREE.get(ref_1letter) if ref_1letter else None
                actual = chain_residue_map[pos]
                if expected:
                    n_checked += 1
                    actual_canonical = NONSTANDARD_TO_CANONICAL.get(actual, actual)
                    if actual_canonical != expected:
                        n_mismatched += 1
                        mismatch_positions.add(pos)
                        mismatch_details.append(
                            f'  pos {pos}: score file expects {ref_1letter} ({expected}), chain has {actual}')

            if n_checked:
                identity_pct = 100 * (n_checked - n_mismatched) / n_checked
                if n_mismatched == 0:
                    validation_log.append(f'Sequence identity: {n_checked}/{n_checked} positions match chain {chain_id} (100%) ✓')
                else:
                    validation_log.append(f'Sequence identity: {n_checked - n_mismatched}/{n_checked} positions match '
                                           f'chain {chain_id} ({identity_pct:.0f}%)')
                    preview = mismatch_details[:10]
                    suffix  = f'\n  ... and {len(mismatch_details) - 10} more' if len(mismatch_details) > 10 else ''
                    validation_log.append('\n'.join(preview) + suffix)
                    if identity_pct < 80:
                        serious_issue = True
                        validation_log.append(f'  *** LOW SEQUENCE IDENTITY — more than 20% of positions do not match chain {chain_id}.'
                                              f' You may have selected the wrong chain or gene. ***')

            # --- Residue number offset detection ---
            # Offset dialog runs now (affects coloring); result appended to deferred log.
            if n_checked == 0 or identity_pct < 90:
                print('Scanning for residue number offset (e.g. signal peptide, PDB renumbering)...')
                best_offset, best_matches, n_checkable = find_best_offset(ref_aa_by_pos, chain_residue_map)
                current_matches = (n_checked - n_mismatched) if n_checked else 0
                improvement = best_matches - current_matches

                if best_offset != 0 and n_checkable > 0 and improvement > max(5, 0.05 * n_checkable):
                    best_pct = 100 * best_matches / n_checkable
                    choice, ok = QInputDialog.getItem(
                        parent, 'Residue Number Offset Detected',
                        f'Applying an offset of {best_offset:+d} to score-file positions\n'
                        f'improves sequence identity from {identity_pct:.0f}% → {best_pct:.0f}%.\n\n'
                        f'This commonly happens when the PDB model starts at a different\n'
                        f'residue than the canonical sequence (e.g. signal peptide removed).\n\n'
                        f'Apply offset {best_offset:+d}?',
                        ['Yes — apply offset', 'No — continue without offset'],
                        0, False)
                    if ok and choice.startswith('Yes'):
                        normalized_values = {pos + best_offset: val
                                             for pos, val in normalized_values.items()}
                        ref_aa_by_pos     = {pos + best_offset: aa
                                             for pos, aa in ref_aa_by_pos.items()}
                        # Confirm identity after offset; rebuild mismatch_positions for coloring
                        common_off = set(normalized_values) & set(chain_residue_map)
                        mismatch_positions = set()
                        n_mis_off = 0
                        for pos in common_off:
                            ref_1l = ref_aa_by_pos.get(pos)
                            exp    = AA_ONE_TO_THREE.get(ref_1l) if ref_1l else None
                            if exp and NONSTANDARD_TO_CANONICAL.get(chain_residue_map[pos], chain_residue_map[pos]) != exp:
                                n_mis_off += 1
                                mismatch_positions.add(pos)
                        n_off = len(common_off)
                        id_off = 100 * (n_off - n_mis_off) / n_off if n_off else 0
                        validation_log.append(f'Offset {best_offset:+d} applied: {n_off - n_mis_off}/{n_off} '
                                              f'positions match ({id_off:.0f}%)')
                        if id_off >= 80:
                            serious_issue = False
                    else:
                        validation_log.append(f'Offset {best_offset:+d} suggested but not applied.')
                else:
                    validation_log.append('No residue number offset improvement found.')

            # --- Gate coloring on validation result ---
            if serious_issue:
                proceed, ok = QInputDialog.getItem(
                    parent, 'Validation Failed — Skip Coloring?',
                    f'Serious issues were found for chain {chain_id} (see log above).\n\n'
                    f'Coverage: {coverage_pct:.0f}%   |   '
                    f'Sequence identity: {identity_pct:.0f}%\n\n'
                    f'Coloring with bad data will produce a misleading result.\n'
                    f'Proceed anyway?',
                    ['No — skip this chain', 'Yes — color anyway (not recommended)'],
                    0, False)
                if not ok or proceed.startswith('No'):
                    print(f'Skipping coloring for chain {chain_id}.')
                    skip_chain = True

        if not skip_chain:
            print(f'Coloring chain {chain_id}...')

            if not existing_models:
                run(session, f'show /{chain_id} cartoons')  #Shows cartoon for chain
                run(session, f'hide /{chain_id} & protein atoms')     #Hides atom representation for protein only
                run(session, f'hide /{chain_id} & protein bonds')     #Hides bond/stick representation for protein only
            run(session, f'surface /{chain_id} & protein')  #Show surface so score colors apply to it; user can hide if not wanted
            run(session, f'color /{chain_id} & protein gray target abcs') #Colors protein residues grey first (cartoons, atoms, surface), excludes pseudobonds (e.g. H-bonds)

            #this block does the coloring; mismatched positions are left gray
            for residue, value in normalized_values.items():
                if residue in mismatch_positions:
                    continue
                color = get_color(value) #Gets color from color map
                hex_color = rgb_to_hex(color[0], color[1], color[2])
                run(session, f'color /{chain_id}:{residue} {hex_color} target abcs') #Colors cartoons, atoms, and surface

            # Print deferred validation summary now that coloring output has finished
            print(f'--- Validation summary for chain {chain_id} ---')
            for line in validation_log:
                print(line)

        again, ok = QInputDialog.getText(
            parent, 'Add Another?', 'Color another gene/chain? (y/n):')
        if not ok or again.strip().lower() != 'y':
            break

    run(session, 'lighting flat')
    print('Done!')

main()
