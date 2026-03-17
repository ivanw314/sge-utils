"""Generate a VCF file from an SGE allscores.tsv for submission to Ensembl VEP.

This script converts the SNV rows of an *allscores.tsv file into a minimal
VCF (GRCh38) that can be uploaded to the Ensembl VEP web tool to obtain
pathogenicity predictor scores (AlphaMissense, REVEL, CADD, SpliceAI).
The downloaded VEP Excel output can then be placed in the pipeline input
directory and will be auto-detected as a *{gene}*vep* file.

Usage:
    python make_vcf.py <allscores_file> [output_vcf]

Arguments:
    allscores_file  Path to the *allscores.tsv file (required)
    output_vcf      Path for the output VCF file (optional).
                    If omitted, the VCF is written to the same directory as
                    the input file, replacing the .tsv extension with .vcf.

Example:
    python make_vcf.py ./data/BRCA1/BRCA1.allscores.tsv
    python make_vcf.py ./data/BRCA1/BRCA1.allscores.tsv ./data/BRCA1/BRCA1.vcf

Next steps after running this script
-------------------------------------
1. Go to https://www.ensembl.org/Tools/VEP (ensure theg enome build is set to GRCh38).

2. Upload the generated .vcf file.

3. Under "Identifiers and frequency data", enable:
      - Gene symbol
      - MANE

4. Under "Pathogenicity predictions", enable:
      - AlphaMissense
      - CADD
      - REVEL
    - SpliceAI

5. Run VEP, then download the results as Excel (.xlsx).

6. Rename the downloaded file to match the pattern: {gene}.vep.xlsx
   e.g. BRCA1.vep.xlsx

7. Place it in the same input directory as your allscores.tsv.
   The pipeline will detect and load it automatically.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an SGE allscores.tsv to a VCF for Ensembl VEP."
    )
    parser.add_argument(
        "allscores",
        type=Path,
        help="Path to the *allscores.tsv file",
    )
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=None,
        help="Output VCF path (default: same directory as input, .vcf extension)",
    )
    return parser.parse_args()


def make_vcf(allscores_path: Path, output_path: Path):
    df = pd.read_csv(allscores_path, sep="\t")

    # SNVs only — deletions cannot be submitted to VEP as single-position entries
    df = df.loc[df["ref"].str.len() == 1].copy()

    df = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    df = df.sort_values(["chrom", "pos"])

    with open(output_path, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        f.write("##reference=GRCh38\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in df.iterrows():
            f.write(
                f"{row['chrom']}\t{row['pos']}\t.\t"
                f"{row['ref'].upper()}\t{row['alt'].upper()}\t.\tPASS\t.\n"
            )

    print(f"VCF written to:   {output_path}")
    print(f"Total variants:   {len(df)}")
    print()
    print("Next steps:")
    print("  1. Upload the VCF to https://www.ensembl.org/Tools/VEP (ensure genome build is GRCh38)")
    print("  2. Enable pathogenicity predictions: AlphaMissense, CADD, REVEL, SpliceAI")
    print("  3. Download results as Excel (.xlsx) or plain text (.txt) — both are supported")
    stem = allscores_path.stem.split('allscores')[0].rstrip('.')
    print(f"  4. Rename to {stem}.vep.xlsx  (or {stem}.vep.txt)")
    print("  5. Place in your pipeline input directory")


def main():
    args = parse_args()

    if not args.allscores.is_file():
        sys.exit(f"Error: file not found: {args.allscores}")

    if "allscores" not in args.allscores.name:
        print(f"Warning: '{args.allscores.name}' does not look like an allscores file.")

    output = args.output
    if output is None:
        output = args.allscores.with_suffix(".vcf")

    make_vcf(args.allscores, output)


if __name__ == "__main__":
    main()
