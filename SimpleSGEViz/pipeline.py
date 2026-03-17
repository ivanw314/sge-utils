"""SGE Visualization Pipeline

Usage:
    python pipeline.py <input_dir> <output_dir> [--format html|png|svg] [--excel]

The input directory must contain:
    *allscores.tsv      SNV + deletion fitness scores (combined file)
    *modelparams.tsv    SGE model thresholds
    *snvcounts.tsv      Per-replicate SNV counts
    *delcounts.tsv      Per-replicate deletion counts

Optional (figures generated only if detected):
    *{gene}*gnomAD*     gnomAD allele frequencies (CSV or Excel)
    *{gene}*Regeneron*  Regeneron allele frequencies (CSV or Excel)
    *{gene}*editrates*  Library edit rates (TSV with target_rep + edit_rate columns)
    *{gene}*targets*    Library targets TSV (columns: editstart, editstop; used as library amplicons)
    *{gene}*vep*        VEP Excel output (.xlsx) with AlphaMissense, REVEL, CADD, SpliceAI scores

Outputs (saved to output_dir):
    {gene}_histogram_stripplot    Score distribution histogram + strip plot
    {gene}_correlation_heatmap    Replicate Pearson r heatmap
    {gene}_scores_across_gene     Per-exon scatter plot of scores vs. genomic position
    {gene}_aa_heatmap             Amino acid substitution heatmap (if amino_acid_change column present)
    {gene}_predictor_scatter      Predictor vs. fitness score scatter panels (if predictor columns present)
    {gene}_clinvar_strip          ClinVar classification strip plot (if *{gene}*ClinVar*SNV* file present)
    {gene}_clinvar_roc            ROC curve for SGE score B/LB vs P/LP classification (if ClinVar file present)
    {gene}_maf_vs_score           Allele frequency vs. score heatmap (if AF files present)
    {gene}_edit_rate_barplot      Library edit rate bar plot by target (if *{gene}*editrates* file present)
    {gene}_exon_cartoon           Exon structure cartoon (fetched from Ensembl; no *targets* file)
    {gene}_library_cartoon        Exon + library design cartoon (fetched from Ensembl + *targets* file)
    {gene}_data.xlsx              Multi-sheet Excel workbook (if --excel flag is set)

PNG and SVG output require vl-convert-python (pip install vl-convert-python).
Excel output requires openpyxl (pip install openpyxl).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from sgeviz import io, process
from sgeviz.figures import aa_heatmap, clinvar_strip, correlation, edit_rate_barplot, gene_cartoon, histogram_strip, maf_score, predictor_scatter, scores_gene


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate standard SGE visualization figures from raw score files."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing input TSV files",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory to write output figures",
    )
    parser.add_argument(
        "--format",
        choices=["html", "png", "svg"],
        default="png",
        help="Output format for figures (default: html). "
             "PNG/SVG require vl-convert-python.",
    )
    parser.add_argument(
        "--excel",
        action="store_true",
        default=False,
        help="Also write a multi-sheet Excel workbook ({gene}_annotated_data.xlsx) "
             "containing scores, thresholds, counts, and any optional datasets. "
             "Requires openpyxl.",
    )
    parser.add_argument(
        "--protein-length",
        type=int,
        default=None,
        metavar="N",
        help="Known full protein length (aa). If the data covers fewer residues, "
             "the x-axis will be extended to this length. If omitted, you will be "
             "prompted interactively for each gene.",
    )
    parser.add_argument(
        "--px-per-aa",
        type=int,
        default=4,
        metavar="N",
        help="Pixels per amino acid column in the AA heatmap (default: 4). "
             "Reduce to produce a narrower figure.",
    )
    parser.add_argument(
        "--gene-name",
        type=str,
        default=None,
        metavar="NAME",
        help="Override the gene name used in figure titles and output filenames. "
             "Cannot be used when multiple gene datasets are detected.",
    )
    parser.add_argument(
        "--assembly",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38",
        help="Genome assembly for Ensembl coordinate fetching (default: GRCh38).",
    )
    parser.add_argument(
        "--exon-color",
        type=str,
        default=None,
        metavar="HEX",
        help="Override exon color in gene cartoons (e.g. '#2E86C1'). "
             "Defaults to the value in the cartoon metadata file, or '#2E86C1'.",
    )
    parser.add_argument(
        "--lib-color",
        type=str,
        default=None,
        metavar="HEX",
        help="Override library amplicon color in gene cartoons (e.g. '#888888'). "
             "Defaults to the value in the cartoon metadata file, or '#888888'.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"Error: input directory not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fmt = args.format

    # --- Discover genes ---
    print(f"Scanning for gene datasets in: {args.input_dir}")
    genes = io.find_genes(args.input_dir)
    print(f"  Found {len(genes)} gene(s): {', '.join(genes)}")

    if args.gene_name is not None:
        if len(genes) > 1:
            sys.exit(
                f"Error: --gene-name cannot be used when multiple gene datasets are detected "
                f"({', '.join(genes)})."
            )
        original = next(iter(genes))
        genes = {args.gene_name: genes[original]}
        print(f"  Gene name overridden: '{original}' -> '{args.gene_name}'")

    # --- Process each gene ---
    for gene, files in genes.items():
        print(f"\n[{gene}] Loading data...")
        scores_df, thresholds = process.load_scores(files)
        scores_df = process.load_vep(files, scores_df)
        counts_df = io.load_counts(files)
        print(f"  {len(scores_df)} variants loaded")

        # --- Fetch exon coords from Ensembl (aa_exon_df needed for heatmap) ---
        cartoon_data = io.load_cartoon(files)
        if cartoon_data is None:
            print(f"[{gene}] Querying Ensembl for canonical transcript...")
            try:
                tx_info = io.get_canonical_transcript(gene, assembly=args.assembly)
                canonical_flag = " [canonical]" if tx_info["is_canonical"] else " [longest coding]"
                print(
                    f"  Auto-selected: {tx_info['transcript_id']}"
                    f"  ({tx_info['biotype']}, {tx_info['n_exons']} exons, "
                    f"{tx_info['strand']}-strand){canonical_flag}"
                )
                raw = input(
                    "  Press Enter to use this transcript, "
                    "or enter a different Ensembl transcript ID: "
                ).strip()
                chosen_tx = raw if raw else None
                cartoon_data = io.fetch_exon_coords(
                    gene,
                    transcript_id=chosen_tx,
                    assembly=args.assembly,
                    _raw_data=tx_info["_raw_data"],
                )
            except (ValueError, ConnectionError) as exc:
                print(f"[{gene}] Could not fetch exon coords: {exc}")

        aa_exon_df = None
        inferred_protein_length = None
        if cartoon_data is not None:
            _exon_df, _, _meta_df = cartoon_data
            try:
                aa_exon_df, inferred_protein_length = io.exon_genomic_to_aa(_exon_df, _meta_df)
            except Exception as exc:
                print(f"[{gene}] Could not convert exon coords to AA positions: {exc}")

        # --- Resolve protein length ---
        # Priority: --protein-length flag > inferred from Ensembl CDS > interactive prompt
        if args.protein_length is not None:
            protein_length = args.protein_length
        elif inferred_protein_length is not None:
            protein_length = inferred_protein_length
            print(f"  Protein length inferred from Ensembl CDS: {protein_length} aa")
        else:
            raw = input(f"  Protein length for {gene} (press Enter to estimate from data): ").strip()
            protein_length = int(raw) if raw else None

        print(f"[{gene}] Generating figures (format: {fmt})")

        hist, strip = histogram_strip.make_figures(scores_df, thresholds, gene=gene)
        io.save_figure(
            histogram_strip.combine(hist, strip),
            args.output_dir / f"{gene}_histogram_stripplot.{fmt}",
        )

        r_df = correlation.compute_correlations(counts_df)
        io.save_figure(
            correlation.make_heatmap(r_df, gene=gene),
            args.output_dir / f"{gene}_correlation_heatmap.{fmt}",
        )

        io.save_figure(
            scores_gene.make_plot(scores_df, thresholds, gene=gene),
            args.output_dir / f"{gene}_scores_across_gene.{fmt}",
        )

        if "amino_acid_change" in scores_df.columns:
            domains_path = files.get("domains")
            if domains_path is not None:
                print(f"[{gene}] Domain file detected: {domains_path.name}")
            else:
                print(f"[{gene}] No domain file detected (place a *{gene}*domain* file in input dir to enable).")
            io.save_figure(
                aa_heatmap.make_plot(
                    scores_df, gene=gene, thresholds=thresholds,
                    domains_path=domains_path,
                    protein_length=protein_length,
                    px_per_aa=args.px_per_aa,
                    aa_exon_df=aa_exon_df,
                ),
                args.output_dir / f"{gene}_aa_heatmap.{fmt}",
            )
        else:
            print(f"[{gene}] No amino_acid_change column, skipping AA heatmap.")

        pred_plot = predictor_scatter.make_plot(scores_df, thresholds, gene=gene)
        if pred_plot is not None:
            io.save_figure(
                pred_plot,
                args.output_dir / f"{gene}_predictor_scatter.{fmt}",
            )
        else:
            print(f"[{gene}] No predictor score columns found, skipping predictor scatter.")

        clinvar_df = process.load_clinvar(files, scores_df)
        if clinvar_df is not None:
            io.save_figure(
                clinvar_strip.make_strip(clinvar_df, thresholds, gene=gene),
                args.output_dir / f"{gene}_clinvar_strip.{fmt}",
            )
            roc = clinvar_strip.make_roc(clinvar_df, gene=gene)
            if roc is not None:
                io.save_figure(
                    roc,
                    args.output_dir / f"{gene}_clinvar_roc.{fmt}",
                )
            else:
                print(f"[{gene}] Insufficient B/LB and P/LP variants for ROC curve.")
        else:
            print(f"[{gene}] No ClinVar file found, skipping ClinVar figures.")

        maf_df = process.load_allele_freqs(files, scores_df)
        if maf_df is not None:
            io.save_figure(
                maf_score.make_plot(maf_df, gene=gene),
                args.output_dir / f"{gene}_maf_vs_score.{fmt}",
            )
        else:
            print(f"[{gene}] No allele frequency files found, skipping MAF figure.")

        targets_lib_df = io.load_targets(files)
        if targets_lib_df is not None:
            print(f"[{gene}] Targets file detected: {files['targets'].name} ({len(targets_lib_df)} amplicons)")
        if cartoon_data is not None:
            exon_df, lib_df_cartoon, meta_df = cartoon_data
            lib_df = targets_lib_df if targets_lib_df is not None else lib_df_cartoon
            if lib_df is not None and not lib_df.empty:
                cartoon_chart = gene_cartoon.make_library_cartoon(
                    exon_df, lib_df, meta_df,
                    exon_color=args.exon_color,
                    lib_color=args.lib_color,
                )
                cartoon_name = f"{gene}_library_cartoon"
            else:
                cartoon_chart = gene_cartoon.make_exon_cartoon(
                    exon_df, meta_df,
                    exon_color=args.exon_color,
                )
                cartoon_name = f"{gene}_exon_cartoon"
            io.save_figure(
                cartoon_chart,
                args.output_dir / f"{cartoon_name}.{fmt}",
            )
        else:
            print(f"[{gene}] No exon coordinates available, skipping gene cartoon.")

        edit_rates_df = io.load_edit_rates(files)
        if edit_rates_df is not None:
            io.save_figure(
                edit_rate_barplot.make_plot(edit_rates_df, gene=gene),
                args.output_dir / f"{gene}_edit_rate_barplot.{fmt}",
            )
        else:
            print(f"[{gene}] No edit rates file found, skipping edit rate bar plot.")

        if args.excel:
            print(f"[{gene}] Writing Excel workbook...")
            thresh_df = pd.DataFrame({
                "non_functional_threshold": [thresholds[0]],
                "functional_threshold": [thresholds[1]],
            })

            # Build merged scores sheet: left-join optional annotations by pos_id
            merged_scores = scores_df.copy()
            if clinvar_df is not None:
                merged_scores = pd.merge(
                    merged_scores,
                    clinvar_df[["pos_id", "Germline classification"]],
                    on="pos_id",
                    how="left",
                )
            if maf_df is not None:
                af_wide = (
                    maf_df[["pos_id", "dataset", "Allele Frequency"]]
                    .pivot_table(
                        index="pos_id",
                        columns="dataset",
                        values="Allele Frequency",
                        aggfunc="first",
                    )
                    .reset_index()
                )
                af_wide.columns.name = None
                af_wide = af_wide.rename(columns={
                    "gnomAD": "gnomad_af",
                    "Regeneron": "regeneron_af",
                })
                merged_scores = pd.merge(merged_scores, af_wide, on="pos_id", how="left")

            sheets = {"scores": merged_scores, "thresholds": thresh_df, "counts": counts_df}
            if edit_rates_df is not None:
                sheets["edit_rates"] = edit_rates_df
            io.save_excel(sheets, args.output_dir / f"{gene}_data.xlsx")

    print(f"\nDone. Figures saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
