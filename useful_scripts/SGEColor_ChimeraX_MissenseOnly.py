"""
SGEColor_ChimeraX_MissenseOnly.py
==================================
Colors a protein ribbon structure in ChimeraX by per-residue SGE (Saturating Genome Editing)
scores, using only missense variants. Scores are aggregated per amino acid position and mapped
onto a white (benign/neutral) → red (damaging) color scale clamped to [-0.2, 0].

USAGE
-----
In the ChimeraX command line:
    runscript /path/to/SGEColor_ChimeraX_MissenseOnly.py

INTERACTIVE DIALOGS (shown at runtime)
---------------------------------------
1. PDB ID          — prompted only if no structure is currently loaded.
                     If a model is already open, that structure is used as-is.
2. SGE score file  — file picker for an Excel (.xlsx) SGE score table.
                     Supported formats: .xlsx (must contain a sheet named 'scores'),
                     .tsv (tab-separated), .csv (comma-separated).
                     Required columns: variant_qc_flag, consequence,
                       amino_acid_change, score
                     RNA_score (optional column) used by RNA filter below
3. Chain selection — dropdown populated from the chains in the loaded structure.
4. RNA filter      — (optional) enter a numeric threshold to exclude variants
                     where RNA_score is below that value. Cancel to skip.
5. Add another?    — repeat steps 2–4 to color additional chains in one run.

SURFACE COLORING
-----------------
The script colors atoms, cartoon, and surface (target abcs). For surface coloring
to apply, the surface must already be visible before running the script. To show
the surface first, run in the ChimeraX command line or select "Show" under "Surfaces" in the "Molecule Display" tab:
    surface

CONFIGURATION (edit at top of script)
---------------------------------------
  analysis_type  'med' | 'mean' | 'min'   Aggregation method per residue (default: 'med')
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

analysis_type = 'med' #Specify whether to color based on 'med', 'mean', or 'min' SGE scores
show_legend = False #Whether to show the legend figure
save_legend = False #Whether to save the legend figure
dna_style = 'stubs' #ssDNA display style: 'stubs', 'slab', 'fill', or 'atoms'. ('ladder' requires dsDNA and won't work here)


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

    df = df.loc[df['variant_qc_flag'] != 'WARN'] #Filters out variants with WARN flag
    df = df.rename(columns={'consequence': 'Consequence', 'amino_acid_change': 'AAsub', 'score': 'snv_score'})
    df = df.loc[df['Consequence'].str.contains('missense_variant')] #Filters only for missense variants

    if rna_score_threshold is not None: #Optionally filter out variants below the RNA score threshold; NaN rows are kept
        df = df.loc[df['RNA_score'].isna() | (df['RNA_score'] >= rna_score_threshold)]

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
    # First clamp all values between 0 and 1
    clamped_values = {k: min(max(v, -0.2), 0) for k, v in values.items()}
    clamped_values = {k: v for k, v in clamped_values.items() if not pd.isna(v)} #Filters out NA values

    # Get min and max of clamped values
    min_val = min(clamped_values.values())
    max_val = max(clamped_values.values())

    # Avoid division by zero if all values are the same
    if max_val == min_val:
        return {k: 0 for k in values.keys()}

    # Perform normalization
    return {k: (v - min_val) / (max_val - min_val) for k, v in clamped_values.items()}


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
    # Create figure with horizontal proportions
    fig, ax = plt.subplots(figsize=(1, 0.5))
    fig.subplots_adjust(bottom=0.5)

    # Reverse the colormap to match your get_color inversion
    reversed_cmap = custom_cmap.reversed()

    # Now use normal ordering since we reversed the colormap
    norm = plt.Normalize(vmin=-0.2, vmax=0)
    sm = cm.ScalarMappable(cmap=reversed_cmap, norm=norm)
    sm.set_array([])

    # Create horizontal colorbar
    cbar = plt.colorbar(sm, cax=ax, orientation='horizontal')

    # Labels now directly correspond to your data range
    cbar.set_ticks([-0.2, 0])
    cbar.set_label('Score', labelpad=10)

    if show_legend:
        plt.show()
    return fig


def get_color(value): #Gets color for each residue from mean score
    return custom_cmap(1 - value) #Inverts the colors, White is high score and Red is low


def rgb_to_hex(r, g, b): #Converts float RGB (0-1) to hex string for ChimeraX
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


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

    legend = create_colorbar_legend()
    if save_legend:
        save_path, _ = QFileDialog.getSaveFileName(parent, 'Save Legend', '', 'PNG Files (*.png)')
        if save_path:
            legend.savefig(save_path, dpi=500)

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

        # Warn about scored positions absent from the structure
        missing = sorted(set(normalized_values) - set(chain_residue_map))
        if missing:
            print(f'Warning: {len(missing)} scored position(s) not found in chain {chain_id}: {missing}')

        # Warn about positions where the reference amino acid doesn't match the structure
        mismatches = []
        for pos in sorted(set(normalized_values) & set(chain_residue_map)):
            ref_1letter = ref_aa_by_pos.get(pos)
            expected = AA_ONE_TO_THREE.get(ref_1letter) if ref_1letter else None
            actual = chain_residue_map[pos]
            if expected and actual != expected:
                mismatches.append(f'  pos {pos}: score file says {ref_1letter} ({expected}), structure has {actual}')
        if mismatches:
            print(f'Warning: {len(mismatches)} amino acid mismatch(es) in chain {chain_id}:')
            for msg in mismatches:
                print(msg)

        print('Applying colors in ChimeraX...')
        print(f'Coloring chain {chain_id} based on {analysis_type} SGE scores...')

        if not existing_models:
            run(session, f'show /{chain_id} cartoons')  #Shows cartoon for chain
            run(session, f'hide /{chain_id} & protein atoms')     #Hides atom representation for protein only
            run(session, f'hide /{chain_id} & protein bonds')     #Hides bond/stick representation for protein only
        run(session, f'color /{chain_id} & protein gray target abcs') #Colors protein residues grey first (cartoons, atoms, surface), excludes pseudobonds (e.g. H-bonds)

        #this block does the coloring
        for residue, value in normalized_values.items():
            if value == 1:
                hex_color = '#ffffff'
            else:
                color = get_color(value) #Gets color from color map
                hex_color = rgb_to_hex(color[0], color[1], color[2])
            run(session, f'color /{chain_id}:{residue} {hex_color} target abcs') #Colors cartoons, atoms, and surface

    run(session, 'lighting flat')
    print('Done!')

main()
