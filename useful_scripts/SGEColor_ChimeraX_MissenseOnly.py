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
Supported file types: .xlsx (sheet must be named 'scores'), .tsv, .csv

Required columns:
  consequence        str    Only rows containing 'missense_variant' are kept.
  amino_acid_change  str    One-letter ref AA + integer position + one-letter alt AA.
                            Examples: 'A123G' (Ala→Gly at position 123),
                                      'R45W'  (Arg→Trp at position 45).
                            The ref AA (first character) is used to verify the
                            sequence matches the selected PDB chain.
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
1. PDB ID          — prompted only if no structure is currently loaded.
                     If a model is already open, that structure is used as-is.
2. Score settings  — choose aggregation method (median / mean / min), then
                     optionally override score column, color range, and direction.
3. Legend          — choose whether to show and/or save the colorbar legend.
4. SGE score file  — file picker for the score table (see INPUT FILE FORMAT above).
5. Chain selection — dropdown populated from the chains in the loaded structure.
6. RNA filter      — (optional) enter a numeric threshold to exclude variants
                     where RNA_score is below that value. Cancel to skip.
7. Add another?    — repeat steps 4–6 to color additional chains in one run.

SURFACE COLORING
-----------------
The script colors atoms, cartoon, and surface (target abcs). For surface coloring
to apply, the surface must already be visible before running the script. To show
the surface first, run in the ChimeraX command line or select "Show" under "Surfaces" in the "Molecule Display" tab:
    surface

CONFIGURATION (edit at top of script)
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

analysis_type = 'med'   # 'med', 'mean', or 'min' — aggregation method per residue
score_column  = 'score' # column in the score file to use for coloring
clamp_min     = -0.2    # lower bound for color normalization range
clamp_max     = 0.0     # upper bound for color normalization range
high_is_red   = False   # False → high score = white, low = red (default); True → high = red, low = white
show_legend = False     # whether to show the legend figure
save_legend = False     # whether to save the legend figure
dna_style = 'stubs'     # ssDNA display style: 'stubs', 'slab', 'fill', or 'atoms'


def read_scores(file, rna_score_threshold=None): #Reads score file
    ext = os.path.splitext(file)[1].lower()
    if ext == '.xlsx':
        df = pd.read_excel(file, sheet_name='scores')
    elif ext == '.tsv':
        df = pd.read_csv(file, sep='\t')
    elif ext == '.csv':
        df = pd.read_csv(file)
    else:
        raise ValueError(f'Unsupported file type: {ext}. Use .xlsx, .tsv, or .csv')

    if 'variant_qc_flag' in df.columns:
        df = df.loc[df['variant_qc_flag'] != 'WARN'] #Filters out variants with WARN flag
    if score_column not in df.columns:
        raise ValueError(f"score_column '{score_column}' not found. Available columns: {df.columns.tolist()}")
    df = df.rename(columns={'consequence': 'Consequence', 'amino_acid_change': 'AAsub', score_column: 'snv_score'})
    df = df.loc[df['Consequence'].str.contains('missense_variant')] #Filters only for missense variants

    if rna_score_threshold is not None: #Optionally filter out variants below the RNA score threshold; NaN rows are kept

        columns = df.columns.tolist()
        if 'RNA_score' in columns:
            df = df.loc[df['RNA_score'].isna() | (df['RNA_score'] >= rna_score_threshold)]
        elif 'RNAscore' in columns:
            df = df.loc[df['RNAscore'].isna() | (df['RNAscore'] >= rna_score_threshold)]
        else:
            raise ValueError(
                "rna_score_threshold was specified but neither 'RNA_score' nor 'RNAscore' "
                "column was found in the score file. Available columns: "
                + str(columns)
            )

    df['AApos'] = df['AAsub'].transform(lambda x: int(x[1:-1])) #Creates new amino acid position column

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
                      if chain_residue_map.get(pos + offset) == exp)
        if matches > best_matches:
            best_matches, best_offset = matches, offset
    return best_offset, best_matches, len(scored)


def get_score_config(session):
    """Ask aggregation method and legend preference (always); optionally override score_column / clamp range / color direction."""
    global analysis_type, score_column, clamp_min, clamp_max, high_is_red, show_legend, save_legend
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
        # Score column
        col, ok = QInputDialog.getText(
            parent, 'Score Column',
            'Score column name in the input file:',
            text=score_column)
        if ok and col.strip():
            score_column = col.strip()

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


def get_gene_configs(session, available_chains):
    parent = session.ui.main_window
    gene_configs = []
    while True:
        file_path, _ = QFileDialog.getOpenFileName(
            parent, 'Select SGE Score File', '',
            'Score Files (*.xlsx *.tsv *.csv);;Excel (*.xlsx);;TSV (*.tsv);;CSV (*.csv)')
        if not file_path:
            break
        chain_id, ok = QInputDialog.getItem(
            parent, 'Select Chain',
            f'Chain for {os.path.basename(file_path)}:',
            available_chains, 0, False)
        if not ok:
            break
        rna_threshold, ok = QInputDialog.getDouble(
            parent, 'RNA Score Filter (optional)',
            'Exclude variants where RNA_score is below (cancel to skip):',
            0.0, -10.0, 10.0, 3)
        rna_score_threshold = rna_threshold if ok else None
        gene_configs.append((file_path, chain_id, rna_score_threshold))
        again, ok = QInputDialog.getText(
            parent, 'Add Another?', 'Color another gene/chain? (y/n):')
        if not ok or again.strip().lower() != 'y':
            break
    if not gene_configs:
        raise ValueError('No data files specified — script aborted.')
    return gene_configs


def main():  # 'session' is injected as a global by ChimeraX at runtime via runscript
    parent = session.ui.main_window

    get_score_config(session)  # optionally override score_column / clamp range / color direction

    legend = create_colorbar_legend()
    if save_legend:
        save_path, _ = QFileDialog.getSaveFileName(parent, 'Save Legend', '', 'PNG Files (*.png)')
        if save_path:
            legend.savefig(save_path, dpi=500)
    plt.close(legend)

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

    gene_configs = get_gene_configs(session, available_chains)

    for file_path, chain_id, rna_score_threshold in gene_configs:
        label = os.path.basename(file_path)
        print(f'Reading SGE scores from {label}...')
        raw_scores = read_scores(file_path, rna_score_threshold)
        ref_aa_by_pos = raw_scores.drop_duplicates('AApos').set_index('AApos')['AAsub'].apply(lambda x: x[0]).to_dict()

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
        serious_issue = False
        n_scored  = len(normalized_values)
        common    = set(normalized_values) & set(chain_residue_map)
        missing   = sorted(set(normalized_values) - set(chain_residue_map))

        # Coverage
        coverage_pct = 100 * len(common) / n_scored if n_scored else 0
        print(f'Coverage: {len(common)}/{n_scored} scored positions found in chain {chain_id} ({coverage_pct:.0f}%)')
        if common:
            range_ok, ok = QInputDialog.getItem(
                parent, 'Check Residue Range',
                f'Scored positions matched chain {chain_id} residues {min(common)}–{max(common)}.\n'
                f'Does this range look correct?',
                ['Yes — continue', 'No — skip this chain'],
                0, False)
            if ok and range_ok.startswith('No'):
                print(f'Skipping coloring for chain {chain_id} — residue range rejected.')
                continue
        if missing:
            preview = missing[:10]
            suffix  = f' ... and {len(missing) - 10} more' if len(missing) > 10 else ''
            print(f'  Note: {len(missing)} scored position(s) not found in chain {chain_id}: {preview}{suffix}')

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
                if actual != expected:
                    n_mismatched += 1
                    mismatch_positions.add(pos)
                    mismatch_details.append(
                        f'  pos {pos}: score file expects {ref_1letter} ({expected}), chain has {actual}')

        if n_checked:
            identity_pct = 100 * (n_checked - n_mismatched) / n_checked
            if n_mismatched == 0:
                print(f'Sequence identity: {n_checked}/{n_checked} positions match chain {chain_id} (100%) ✓')
            else:
                print(f'Sequence identity: {n_checked - n_mismatched}/{n_checked} positions match '
                      f'chain {chain_id} ({identity_pct:.0f}%)')
                preview = mismatch_details[:10]
                suffix  = f'\n  ... and {len(mismatch_details) - 10} more' if len(mismatch_details) > 10 else ''
                print('\n'.join(preview) + suffix)
                if identity_pct < 80:
                    serious_issue = True
                    print(f'  *** LOW SEQUENCE IDENTITY — more than 20% of positions do not match chain {chain_id}.'
                          f' You may have selected the wrong chain or gene. ***')

        # --- Residue number offset detection ---
        # Trigger whenever identity is imperfect; also catches the case where coverage is low
        # because positions simply don't exist at offset 0.
        if n_checked == 0 or identity_pct < 90:
            print('Scanning for residue number offset (e.g. signal peptide, PDB renumbering)...')
            best_offset, best_matches, n_checkable = find_best_offset(ref_aa_by_pos, chain_residue_map)
            current_matches = (n_checked - n_mismatched) if n_checked else 0
            improvement = best_matches - current_matches

            if best_offset != 0 and n_checkable > 0 and improvement > max(5, 0.05 * n_checkable):
                best_pct = 100 * best_matches / n_checkable
                print(f'  Possible offset detected: {best_offset:+d} '
                      f'(identity {identity_pct:.0f}% → {best_pct:.0f}%)')
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
                        if exp and chain_residue_map[pos] != exp:
                            n_mis_off += 1
                            mismatch_positions.add(pos)
                    n_off = len(common_off)
                    id_off = 100 * (n_off - n_mis_off) / n_off if n_off else 0
                    print(f'  After offset {best_offset:+d}: {n_off - n_mis_off}/{n_off} '
                          f'positions match ({id_off:.0f}%)')
                    # Clear serious_issue if identity is now acceptable
                    if id_off >= 80:
                        serious_issue = False
            else:
                print('  No offset improvement found.')

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
                continue

        print('Applying colors in ChimeraX...')
        print(f'Coloring chain {chain_id} based on {analysis_type} SGE scores...')

        if not existing_models:
            run(session, f'show /{chain_id} cartoons')  #Shows cartoon for chain
            run(session, f'hide /{chain_id} & protein atoms')     #Hides atom representation for protein only
            run(session, f'hide /{chain_id} & protein bonds')     #Hides bond/stick representation for protein only
        run(session, f'color /{chain_id} & protein gray target abcs') #Colors protein residues grey first (cartoons, atoms, surface), excludes pseudobonds (e.g. H-bonds)

        #this block does the coloring; mismatched positions are left gray
        for residue, value in normalized_values.items():
            if residue in mismatch_positions:
                continue
            if value == 1:
                hex_color = '#ffffff'
            else:
                color = get_color(value) #Gets color from color map
                hex_color = rgb_to_hex(color[0], color[1], color[2])
            run(session, f'color /{chain_id}:{residue} {hex_color} target abcs') #Colors cartoons, atoms, and surface

    run(session, 'lighting flat')
    print('Done!')

main()
