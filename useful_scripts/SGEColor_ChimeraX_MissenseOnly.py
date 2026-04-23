from chimerax.core.commands import run
import subprocess
import sys

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

#This script generates colored ribbon structures in ChimeraX based on SGE scores.
#ChimeraX port of BCDX2_colorPyMOL_MIS_only.py

#To run: open ChimeraX, then use the command:
#   runscript /path/to/BCDX2_colorChimeraX_MIS_only.py


analysis_type = 'med' #Specify whether to color based on 'med' or 'min' SGE scores

rad51d_version_date = '20260407' #Version date for RAD51D SGE score file, should match version date used in BCDX2_MakeFinalDataTable.ipynb
xrcc2_version_date = '20260122' #Version date for XRCC2 SGE score file, should match version date used in BCDX2_MakeFinalDataTable.ipynb

sge_data = {'RAD51D': f'/Users/ivan/Documents/GitHub/RAD51D_figs/Data/final_tables/supplementary_file_1_RAD51D_SGE_final_table_{rad51d_version_date}.xlsx',
     'XRCC2': f'/Users/ivan/Documents/GitHub/RAD51D_figs/Data/final_tables/supplementary_file_1_XRCC2_SGE_final_table_{xrcc2_version_date}.xlsx'
    }

pdb = '8GBJ' #PDB structure to use

show_legend = False #Whether to show the legend figure
save_legend = False #Whether to save the legend figure
filter_rna_low = True #Whether to filter out RAD51D variants with RNA_consequence == 'low'
dna_style = 'stubs' #ssDNA display style: 'stubs', 'slab', 'fill', or 'atoms'. ('ladder' requires dsDNA and won't work here)

#The chains are specified as follows based on user-provided PDB
if pdb == '8OUZ':
    chains = {'RAD51D': 'C',
            'XRCC2': 'D'}
elif pdb == '8GBJ':
    chains = {'RAD51D': 'D',
            'XRCC2': 'X'}


def region_residues(gene): #Takes gene input and creates the respective residue numbers
    region_residues = [] #Initializes empty list for region residues

    if gene == 'RAD51D':
        for i in range(1, 329):
            region_residues.append(i)
    elif gene == 'XRCC2':
        for i in range(1, 281):
            region_residues.append(i)

    return region_residues

def read_scores(file, region_resi, gene): #Reads score file
    df = pd.read_excel(file, sheet_name='scores') #Reads in score file

    df = df.loc[df['variant_qc_flag'] != 'WARN'] #Filters out variants with WARN flag
    df = df.rename(columns = {'consequence': 'Consequence', 'amino_acid_change': 'AAsub', 'score': 'snv_score'})
    df = df.loc[df['Consequence'].str.contains('missense_variant')] #Filters only for missense variants

    if gene == 'RAD51D' and filter_rna_low: #Filters out RAD51D variants with low RNA consequence
        df = df.loc[df['RNA_consequence'] != 'low']

    df['AApos'] = df['AAsub'].transform(lambda x: int(x[1:-1])) #Creates new amino acid position column

    df = df.loc[df['AApos'].isin(region_resi)] #Filters for variants in residues in region of interest

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


def main():  # 'session' is injected as a global by ChimeraX at runtime via runscript
    legend = create_colorbar_legend()

    if save_legend:
        legend.savefig('/Users/ivan/Desktop/RAD51D_XRCC2_figs/chimerax_ribbon_legend.png', dpi = 500)

    existing_models = session.models.list()
    if existing_models:
        print(f'Existing session detected ({len(existing_models)} model(s) open). Skipping structure load — applying colors only.')
    else:
        print(f'Loading structure {pdb}...')
        run(session, f'open {pdb} from pdb') #Fetches PDB structure from RCSB

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

    for gene in sge_data.keys():
        chain = chains[gene] #Gets chain for specified gene

        print(f'Reading SGE scores for {gene}...')
        region_resi = region_residues(gene) #Gets region residues
        sge_scores = sge_data[gene] #Gets SGE score file path
        raw_scores = read_scores(sge_scores, region_resi, gene) #Gets raw SGE scores

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

        print('Applying colors in ChimeraX...')
        print(f'Coloring chain {chain} based on {analysis_type} SGE scores...')

        if not existing_models:
            run(session, f'show /{chain} cartoons')  #Shows cartoon for chain
            run(session, f'hide /{chain} & protein atoms')     #Hides atom representation for protein only
            run(session, f'hide /{chain} & protein bonds')     #Hides bond/stick representation for protein only
        run(session, f'color /{chain} & protein gray target abc') #Colors protein residues grey first (cartoons and atoms), excludes pseudobonds (e.g. H-bonds)

        #this block does the coloring
        for residue, value in normalized_values.items():
            if value == 1:
                hex_color = '#ffffff'
            else:
                color = get_color(value) #Gets color from color map
                hex_color = rgb_to_hex(color[0], color[1], color[2])
            run(session, f'color /{chain}:{residue} {hex_color}') #Colors cartoons and atoms

    print('Done!')
main()
