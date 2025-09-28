#!/usr/bin/env python3
import os
import sys
import gzip
import pandas as pd
from pathlib import Path

# -------- Nextflow template variables (expanded at runtime) --------
# Required minimal inputs
PREFIX              = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
BARCODES_PATH       = "${barcodes}"

# Optional tool outputs (empty string if not provided)
VIREO_DONOR_IDS     = "${vireo_donor_ids}"
VIREO_SUMMARY       = "${vireo_summary}"            # not strictly needed; kept for future
DEMUXLET_RESULT     = "${demuxlet_result}"
FREEMUXLET_RESULT   = "${freemuxlet_result}"
SOUPORCELL_TSV      = "${souporcell_tsv}"

# Optional flags + inputs for .h5ad (RNA)
GENERATE_ANNDATA    = "${generate_anndata}"
RNA_MATRIX_10X      = "${rna_matrix}"               # optional: 10x mtx dir path if you want .h5ad

SINGLET  = "singlet"
DOUBLET  = "doublet"
NEGATIVE = "negative"

def _exists(p: str) -> bool:
    return isinstance(p, str) and len(p.strip()) > 0 and Path(p).exists()

def _read_table_maybe_gzip(path: str, **kwargs) -> pd.DataFrame:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh, **kwargs)
    # try TSV then CSV
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except Exception:
        return pd.read_csv(path, **kwargs)

# -------------------- Load per-tool outputs --------------------

def load_vireo(donor_ids_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Vireo donor IDs file typically has columns like: cell, donor_id (headered TSV).
    We use donor_id as assignment; classification is singlet if looks like a donor label,
    otherwise (if missing) leave NA.
    """
    df = _read_table_maybe_gzip(donor_ids_path)
    # try to standardize column names
    cols = {c.lower(): c for c in df.columns}
    cell_col  = cols.get("cell") or cols.get("barcode") or list(df.columns)[0]
    donor_col = cols.get("donor_id") or cols.get("donor") or list(df.columns)[1]

    out = df[[cell_col, donor_col]].copy()
    out.columns = ["Barcode", "vireo"]

    assign = out.copy()
    clas   = out.copy()
    clas["vireo"] = clas["vireo"].apply(
        lambda x: SINGLET if pd.notna(x) and str(x).strip() not in {"", DOUBLET, NEGATIVE} else pd.NA
    )
    return assign, clas

def load_demuxlet(best_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Demuxlet *.best usually contains a BEST field like:
      SNG-IND1, DBL-IND1-IND2, or AMB
    Common columns: BARCODE, BEST (TSV).
    """
    df = _read_table_maybe_gzip(best_path)
    cols = {c.lower(): c for c in df.columns}
    barcode_col = cols.get("barcode") or list(df.columns)[0]
    best_col    = cols.get("best")    or list(df.columns)[1]

    tmp = df[[barcode_col, best_col]].copy()
    tmp.columns = ["Barcode", "BEST"]

    def parse_assignment(s: str):
        if pd.isna(s): return pd.NA
        s = str(s)
        if s.startswith("SNG-"):
            return s.replace("SNG-","",1)
        if s.startswith("DBL-"):
            return DOUBLET
        if s == "AMB":
            return NEGATIVE
        return s

    def parse_class(s: str):
        if pd.isna(s): return pd.NA
        s = str(s)
        if s.startswith("SNG-"): return SINGLET
        if s.startswith("DBL-"): return DOUBLET
        if s == "AMB":          return NEGATIVE
        return pd.NA

    assign = tmp[["Barcode"]].copy()
    assign["demuxlet"] = tmp["BEST"].map(parse_assignment)

    clas = tmp[["Barcode"]].copy()
    clas["demuxlet"] = tmp["BEST"].map(parse_class)

    return assign, clas

def load_freemuxlet(samples_gz_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Freemuxlet *.clust1.samples.gz often has columns:
      BARCODE, CLUST, SNG.BEST, DBL.BEST, etc.
    If SNG.BEST is non-NA => assignment is that donor else DOUBLET.
    """
    df = _read_table_maybe_gzip(samples_gz_path)
    cols = {c.lower(): c for c in df.columns}
    barcode_col = cols.get("barcode") or list(df.columns)[0]
    sngbest_col = cols.get("sng.best") or cols.get("sng_best")
    dblbest_col = cols.get("dbl.best") or cols.get("dbl_best")

    out = df[[barcode_col]].copy()
    out.columns = ["Barcode"]

    def get_assignment(row):
        sng = row.get(sngbest_col) if sngbest_col else None
        dbl = row.get(dblbest_col) if dblbest_col else None
        if pd.notna(sng) and str(sng).strip() not in {"", "NA", "NaN"}:
            return sng
        if pd.notna(dbl) and str(dbl).strip() not in {"", "NA", "NaN"}:
            return DOUBLET
        return pd.NA

    out["freemuxlet"] = df.apply(get_assignment, axis=1)
    assign = out.copy()

    clas = out[["Barcode"]].copy()
    clas["freemuxlet"] = out["freemuxlet"].apply(
        lambda x: SINGLET if pd.notna(x) and x not in {DOUBLET, NEGATIVE} else (DOUBLET if x == DOUBLET else pd.NA)
    )
    return assign, clas

def load_souporcell(tsv_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Souporcell clusters.tsv typically has barcode->cluster/assignment.
    We treat any non-empty non-doublet/non-negative as a singlet assignment.
    """
    df = _read_table_maybe_gzip(tsv_path)
    # heuristics: first column barcode, second the call
    if df.shape[1] < 2:
        raise ValueError("Souporcell TSV appears to have <2 columns; please adjust parser.")
    df = df.iloc[:, :2].copy()
    df.columns = ["Barcode", "souporcell"]

    assign = df.copy()
    clas = df.copy()
    clas["souporcell"] = clas["souporcell"].apply(
        lambda x: SINGLET if pd.notna(x) and str(x) not in {DOUBLET, NEGATIVE, ""} else (str(x) if pd.notna(x) and str(x) in {DOUBLET, NEGATIVE} else pd.NA)
    )
    return assign, clas

# -------------------- main --------------------

def main():
    # base frame of all barcodes (from provided barcodes file; accept 1-col CSV/TSV)
    if not _exists(BARCODES_PATH):
        print("ERROR: barcodes path missing/unreadable.", file=sys.stderr)
        sys.exit(1)

    try:
        bcodes = _read_table_maybe_gzip(BARCODES_PATH, header=None)
    except Exception:
        bcodes = pd.read_csv(BARCODES_PATH, header=None)
    bcodes = bcodes.iloc[:, :1].copy()
    bcodes.columns = ["Barcode"]

    assignment = bcodes.copy()
    classification = bcodes.copy()

    # incrementally merge each tool, if present
    if _exists(VIREO_DONOR_IDS):
        a, c = load_vireo(VIREO_DONOR_IDS)
        assignment = assignment.merge(a, on="Barcode", how="left")
        classification = classification.merge(c, on="Barcode", how="left")

    if _exists(DEMUXLET_RESULT):
        a, c = load_demuxlet(DEMUXLET_RESULT)
        assignment = assignment.merge(a, on="Barcode", how="left")
        classification = classification.merge(c, on="Barcode", how="left")

    if _exists(FREEMUXLET_RESULT):
        a, c = load_freemuxlet(FREEMUXLET_RESULT)
        assignment = assignment.merge(a, on="Barcode", how="left")
        classification = classification.merge(c, on="Barcode", how="left")

    if _exists(SOUPORCELL_TSV):
        a, c = load_souporcell(SOUPORCELL_TSV)
        assignment = assignment.merge(a, on="Barcode", how="left")
        classification = classification.merge(c, on="Barcode", how="left")

    # write CSVs
    out_assign = f"{PREFIX}_genetic_summary_assignment.csv"
    out_class  = f"{PREFIX}_genetic_summary_classification.csv"
    assignment.to_csv(out_assign, index=False)
    classification.to_csv(out_class, index=False)

    # optional: .h5ad, if requested and RNA matrix provided
    make_h5ad = (GENERATE_ANNDATA.lower() == "true") if isinstance(GENERATE_ANNDATA, str) else False
    if make_h5ad and _exists(RNA_MATRIX_10X):
        try:
            import scanpy as sc
            adata = sc.read_10x_mtx(RNA_MATRIX_10X)
            adata.obs = adata.obs.join(assignment.set_index("Barcode"), how="left")
            # fill empties with NEGATIVE for present tools
            tool_cols = [c for c in assignment.columns if c != "Barcode"]
            for col in tool_cols:
                if col in adata.obs.columns:
                    adata.obs[col] = adata.obs[col].astype("string").fillna(NEGATIVE)
            adata.write_h5ad(f"{PREFIX}_genetic_summary.h5ad")
        except Exception as e:
            print(f"WARNING: Failed to write h5ad: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
