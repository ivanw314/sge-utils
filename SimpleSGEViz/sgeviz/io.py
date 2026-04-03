import fnmatch
import json
import urllib.error
import urllib.request
from pathlib import Path

import altair as alt
import pandas as pd


def find_genes(input_dir: Path) -> dict:
    """Discover all gene datasets in the input directory.

    Detects genes by finding all *allscores.tsv files and extracting the gene
    name from each filename. Both dot-separated (e.g. CTCF.allscores.tsv) and
    run-together (e.g. 20260129_RAD51Dallscores.tsv) naming conventions are
    supported. Companion files (*modelparams.tsv, *snvcounts.tsv, *delcounts.tsv)
    are then located using the gene name as a search key.

    Returns a dict mapping gene name -> files dict, e.g.:
        {"RAD51D": {"all_scores": Path(...), "snv_counts": Path(...), ...}}
    """
    allscores_files = sorted(input_dir.glob("*allscores.tsv"))
    if not allscores_files:
        raise FileNotFoundError(f"No '*allscores.tsv' files found in {input_dir}")

    def find_one(*patterns):
        for pattern in patterns:
            matches = list(input_dir.glob(pattern))
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                raise ValueError(
                    f"Multiple files match '{pattern}': "
                    + ", ".join(str(m) for m in matches)
                )
        raise FileNotFoundError(
            f"Could not find any of {patterns} in {input_dir}"
        )

    def find_optional(pattern):
        matches = list(input_dir.glob(pattern))
        return matches[0] if len(matches) == 1 else None

    def find_optional_icase(pattern, exclude=None):
        matches = [
            p for p in input_dir.iterdir()
            if not p.name.startswith("~$")
            and fnmatch.fnmatch(p.name.lower(), pattern.lower())
            and (exclude is None or not fnmatch.fnmatch(p.name.lower(), exclude.lower()))
        ]
        return matches[0] if len(matches) == 1 else None

    genes = {}
    for allscores_path in allscores_files:
        # Handles both GENE.allscores.tsv and GENEallscores.tsv (with optional prefix)
        stem_part = allscores_path.stem.split("_")[-1]
        gene = stem_part.removesuffix(".allscores").removesuffix("allscores")
        genes[gene] = {
            "all_scores": allscores_path,
            "model_params": find_one(f"*{gene}modelparams.tsv", f"*{gene}.modelparams.tsv"),
            "snv_counts": find_one(f"*{gene}snvcounts.tsv", f"*{gene}.snvcounts.tsv"),
            "del_counts": find_one(f"*{gene}delcounts.tsv", f"*{gene}.delcounts.tsv"),
            # Optional allele frequency files (CSV or Excel)
            "gnomad": find_optional(f"*{gene}*gnomAD*"),
            "regeneron": find_optional(f"*{gene}*Regeneron*"),
            # Optional ClinVar SNV file (tab-delimited .txt from ClinVar download)
            "clinvar": find_optional_icase(f"*{gene}*clinvar*snv*"),
            # Optional domain annotation file (Excel or CSV with region_name + aa_residues cols)
            # Exclude files that also contain "cartoon" — those are gene-cartoon files, not domain files.
            "domains": find_optional_icase(f"*{gene}*domain*", exclude=f"*cartoon*"),
            # Optional library edit rates file (*editrates*.tsv)
            "edit_rates": find_optional_icase(f"*{gene}*editrates*"),
            # Optional library targets file (TSV with editstart + editstop columns)
            "targets": find_optional_icase(f"*{gene}*targets*"),
            # Optional gene cartoon file (Excel with exon_coords, metadata; lib_coords sheet no longer needed)
            "cartoon": find_optional_icase(f"*{gene}*cartoon*"),
            # Optional VEP output file (Excel .xlsx or tab-delimited .txt) for predictor scores
            "vep": find_optional_icase(f"*{gene}*vep*"),
        }

    return genes


def load_counts(files: dict) -> pd.DataFrame:
    """Load and merge SNV and deletion counts files, standardizing column names."""
    snv_df = pd.read_csv(files["snv_counts"], sep="\t")
    del_df = pd.read_csv(files["del_counts"], sep="\t")

    del_df["var_type"] = "3bp_del"
    snv_df["start"] = snv_df["pos"]
    snv_df["end"] = snv_df["pos"]
    snv_df["var_type"] = "snv"
    snv_df = snv_df.drop("pos", axis=1)

    df = pd.concat([snv_df, del_df], ignore_index=True)

    # Strip gene prefix from target names (e.g. BARD1_X1B -> 1B, RAD51D_X1 -> 1)
    df["target"] = df["target"].str.replace(r"^[A-Z0-9]+_X", "", regex=True)

    df = df.rename(
        columns={
            "D05_R1": "D05 R1",
            "D05_R2": "D05 R2",
            "D05_R3": "D05 R3",
            "D13_R1": "D13 R1",
            "D13_R2": "D13 R2",
            "D13_R3": "D13 R3",
            "D17_R1": "D17 R1",
            "D17_R2": "D17 R2",
            "D17_R3": "D17 R3",
        }
    )
    return df


def load_cartoon(files: dict):
    """Load a gene cartoon Excel file if present.

    Expected sheets: ``exon_coords`` (required), ``metadata`` (required),
    ``lib_coords`` (optional).

    Returns ``(exon_df, lib_df, metadata_df)`` where ``lib_df`` is ``None``
    when the ``lib_coords`` sheet is absent.  Returns ``None`` if no cartoon
    file was detected.
    """
    path = files.get("cartoon")
    if path is None:
        return None
    xl = pd.ExcelFile(path)
    exon_df = xl.parse("exon_coords")
    lib_df = xl.parse("lib_coords") if "lib_coords" in xl.sheet_names else None
    meta_df = xl.parse("metadata")

    # Coordinates must be in standard genomic order (start < end).
    # Warn and correct any rows where start > end so the user knows their file
    # is non-standard.
    for sheet_name, df in [("exon_coords", exon_df)] + (
        [("lib_coords", lib_df)] if lib_df is not None else []
    ):
        mask = df["start"] > df["end"]
        if mask.any():
            import warnings
            warnings.warn(
                f"Cartoon file '{path.name}': {mask.sum()} row(s) in sheet "
                f"'{sheet_name}' have start > end. Coordinates should be in "
                f"standard genomic order (start < end) regardless of strand. "
                f"The affected rows have been automatically corrected by swapping "
                f"start and end.",
                UserWarning,
                stacklevel=2,
            )
            df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values

    return exon_df, lib_df, meta_df


_SPECIES_MAP = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "zebrafish": "danio_rerio",
    "fly": "drosophila_melanogaster",
    "worm": "caenorhabditis_elegans",
    "yeast": "saccharomyces_cerevisiae",
}


def _fetch_gene_data(gene_symbol: str, species: str, assembly: str) -> dict:
    """Fetch raw gene data from the Ensembl REST API (internal helper)."""
    ens_species = _SPECIES_MAP.get(species.lower(), species.lower())
    base_url = (
        "https://grch37.rest.ensembl.org"
        if assembly.upper() == "GRCH37"
        else "https://rest.ensembl.org"
    )
    url = (
        f"{base_url}/lookup/symbol/{ens_species}/{gene_symbol}"
        "?expand=1&content-type=application/json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        raise ValueError(
            f"Ensembl REST API returned HTTP {e.code} for gene '{gene_symbol}'. "
            "Check that the gene symbol and species are correct."
        ) from e
    except urllib.error.URLError as e:
        raise ConnectionError(
            f"Could not reach Ensembl REST API: {e.reason}. "
            "Check your internet connection."
        ) from e


def _select_transcript(data: dict, gene_symbol: str, transcript_id: str | None):
    """Pick a transcript dict from raw Ensembl gene data (internal helper)."""
    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError(f"No transcripts found for gene '{gene_symbol}'.")
    if transcript_id is not None:
        tx = next((t for t in transcripts if t["id"] == transcript_id), None)
        if tx is None:
            available = [t["id"] for t in transcripts]
            raise ValueError(
                f"Transcript '{transcript_id}' not found for '{gene_symbol}'. "
                f"Available: {available}"
            )
        return tx
    canonical = [t for t in transcripts if t.get("is_canonical") == 1]
    return canonical[0] if canonical else max(
        transcripts,
        key=lambda t: (t.get("Translation") is not None, t.get("length", 0)),
    )


def get_canonical_transcript(
    gene_symbol: str,
    species: str = "human",
    assembly: str = "GRCh38",
) -> dict:
    """Return a summary dict for the canonical transcript of a gene.

    Fetches from the Ensembl REST API and returns a dict with keys:

    - ``transcript_id``: Ensembl transcript ID (e.g. ``"ENST00000357654"``)
    - ``biotype``: transcript biotype (e.g. ``"protein_coding"``)
    - ``n_exons``: number of exons
    - ``is_canonical``: whether Ensembl marks this as the canonical transcript
    - ``strand``: ``"plus"`` or ``"minus"``

    Useful for presenting the auto-selected transcript to a user before
    calling :func:`fetch_exon_coords`.
    """
    data = _fetch_gene_data(gene_symbol, species, assembly)
    tx = _select_transcript(data, gene_symbol, transcript_id=None)
    return {
        "transcript_id": tx["id"],
        "biotype": tx.get("biotype", "unknown"),
        "n_exons": len(tx.get("Exon", [])),
        "is_canonical": bool(tx.get("is_canonical") == 1),
        "strand": "minus" if data["strand"] == -1 else "plus",
        "_raw_data": data,  # passed through to avoid a second API call
    }


def fetch_exon_coords(
    gene_symbol: str,
    transcript_id: str | None = None,
    species: str = "human",
    assembly: str = "GRCh38",
    _raw_data: dict | None = None,
) -> tuple[pd.DataFrame, None, pd.DataFrame]:
    """Fetch exon coordinates for a gene from the Ensembl REST API.

    Returns ``(exon_df, None, meta_df)`` matching the ``load_cartoon()`` output
    format, so the result can be passed directly to
    ``gene_cartoon.make_exon_cartoon()``.

    Args:
        gene_symbol: HGNC gene symbol (e.g. ``"BRCA1"``).
        transcript_id: Ensembl transcript ID to use. If ``None``, the canonical
            transcript is selected automatically.
        species: Common name (e.g. ``"human"``) or Ensembl species name (e.g.
            ``"homo_sapiens"``). Defaults to ``"human"``.
        assembly: Genome assembly — ``"GRCh38"`` (default) or ``"GRCh37"``.
    """
    data = _raw_data or _fetch_gene_data(gene_symbol, species, assembly)
    tx = _select_transcript(data, gene_symbol, transcript_id)
    strand = "minus" if data["strand"] == -1 else "plus"

    # Sort in transcription order (5'→3') so X1 is always the first exon.
    # Plus-strand: ascending genomic start; minus-strand: descending genomic start.
    exons = sorted(tx.get("Exon", []), key=lambda e: e["start"], reverse=(strand == "minus"))
    exon_df = pd.DataFrame([
        {"exon": f"X{i + 1}", "start": e["start"], "end": e["end"]}
        for i, e in enumerate(exons)
    ])

    translation = tx.get("Translation")
    if translation is not None:
        atg_pos = translation["end"] if strand == "minus" else translation["start"]
        stop_pos = translation["start"] if strand == "minus" else translation["end"]
    else:
        atg_pos = int(exon_df["start"].min())
        stop_pos = int(exon_df["end"].max())

    meta_df = pd.DataFrame([
        {"type": "strand", "info": strand},
        {"type": "atg",    "info": atg_pos},
        {"type": "stop",   "info": stop_pos},
    ])

    print(
        f"  Fetched {len(exon_df)} exons for {gene_symbol} "
        f"({tx['id']}, {strand}-strand, {assembly})"
    )
    return exon_df, None, meta_df


def exon_genomic_to_aa(
    exon_df: pd.DataFrame,
    meta_df: pd.DataFrame,
) -> pd.DataFrame:
    """Convert genomic exon coordinates to amino acid residue positions.

    Walks exons in transcription order, computes each exon's CDS overlap
    using the ATG and stop positions in ``meta_df``, and accumulates CDS
    bases to derive per-exon amino acid boundaries.

    Returns a DataFrame with columns ``aa_start`` and ``aa_end`` suitable
    for the exon track in ``aa_heatmap.make_plot()``.  Purely UTR exons
    (no CDS overlap) are omitted.
    """
    meta = dict(zip(meta_df["type"].str.lower(), meta_df["info"].astype(str).str.strip()))
    strand = meta.get("strand", "plus").lower()
    atg_pos = int(float(meta["atg"]))
    stop_pos = int(float(meta["stop"]))

    cds_lo = min(atg_pos, stop_pos)
    cds_hi = max(atg_pos, stop_pos)

    exons = exon_df.sort_values("start", ascending=(strand != "minus"))

    aa_exons = []
    cds_bases_so_far = 0
    for exon_num, (_, exon) in enumerate(exons.iterrows(), start=1):
        ov_start = max(int(exon["start"]), cds_lo)
        ov_end = min(int(exon["end"]), cds_hi)
        # Ensembl uses 1-based fully-closed intervals, so bases = end - start + 1
        cds_bases = max(0, ov_end - ov_start + 1)
        if cds_bases > 0:
            aa_exons.append({
                "exon_num": exon_num,
                "aa_start": round(cds_bases_so_far / 3 + 1, 2),
                "aa_end":   round((cds_bases_so_far + cds_bases) / 3 + 1, 2),
            })
            cds_bases_so_far += cds_bases

    # Subtract 3 bases for the stop codon (included in Ensembl Translation coordinates)
    # before converting to amino acid count.
    protein_length = max(0, cds_bases_so_far - 3) // 3
    return pd.DataFrame(aa_exons, columns=["exon_num", "aa_start", "aa_end"]), protein_length


def load_targets(files: dict) -> "pd.DataFrame | None":
    """Load a library targets TSV file if present.

    Expected columns (at minimum): ``editstart``, ``editstop``.
    Returns a DataFrame with ``start`` and ``end`` columns (genomic coordinates)
    matching the ``lib_df`` format expected by ``gene_cartoon.make_library_cartoon()``,
    or ``None`` if no targets file was detected.
    """
    path = files.get("targets")
    if path is None:
        return None
    df = pd.read_csv(path, sep="\t")
    return df[["editstart", "editstop"]].rename(
        columns={"editstart": "start", "editstop": "end"}
    )


def load_edit_rates(files: dict):
    """Load an edit rates TSV file if present.

    Expected columns: target_rep, edit_rate.
    Returns the raw DataFrame or None if no edit rates file was detected.
    """
    path = files.get("edit_rates")
    if path is None:
        return None
    return pd.read_csv(path, sep="\t")


def save_figure(chart: alt.Chart, path: Path):
    """Save an Altair chart. Format is inferred from the file extension.

    HTML is fully self-contained and interactive.
    PNG and SVG require vl-convert-python to be installed.
    """
    chart.save(str(path))
    print(f"  Saved: {path.name}")


def save_excel(sheets: dict, path: Path):
    """Write a multi-sheet Excel workbook. Requires openpyxl.

    Args:
        sheets: Dict mapping sheet name -> DataFrame (insertion order preserved).
        path: Output path (.xlsx).
    """
    with pd.ExcelWriter(str(path), engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"  Saved: {path.name}")
