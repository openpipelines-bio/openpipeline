"""Complete an ERCC spike-in annotation so it follows the GENCODE/Ensembl layout.

ERCC annotations typically only contain one exon per spike-in with gene_id and
transcript_id. Add the matching gene and transcript features, plus gene/transcript
name and type attributes, so that stricter GTF validators (e.g. BD Rhapsody v4)
accept it. Attributes that are already present in the input are kept, unless
--overwrite_types is set, which replaces existing gene_type and transcript_type values.
"""

import argparse

import pandas as pd
from gtfparse import create_missing_features, read_gtf, write_gtf

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", required=True, help="Input ERCC GTF.")
parser.add_argument("--output", required=True, help="Output GTF.")
parser.add_argument("--gene_type", required=True, help="Value for gene_type.")
parser.add_argument(
    "--transcript_type", required=True, help="Value for transcript_type."
)
parser.add_argument(
    "--overwrite_types",
    action="store_true",
    help="Replace existing gene_type and transcript_type values instead of only filling missing ones.",
)
args = parser.parse_args()

GENE_ATTRIBUTES = ["gene_name", "gene_type"]
TRANSCRIPT_ATTRIBUTES = [
    "gene_id",
    *GENE_ATTRIBUTES,
    "transcript_name",
    "transcript_type",
]
CODING_FEATURES = ["CDS", "start_codon", "stop_codon"]
FEATURE_ORDER = {"gene": 0, "transcript": 1}

df = read_gtf(args.input, result_type="pandas")

# gtfparse uses '' for attributes that are missing on a row; use NA instead
df = df.reindex(columns=df.columns.union(TRANSCRIPT_ATTRIBUTES, sort=False))
df[TRANSCRIPT_ATTRIBUTES] = (
    df[TRANSCRIPT_ATTRIBUTES].astype("string").replace("", pd.NA)
)

is_transcript_level = df["feature"] != "gene"
df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
df["transcript_name"] = df["transcript_name"].fillna(
    df["gene_name"].where(is_transcript_level)
)
if args.overwrite_types:
    df["gene_type"] = args.gene_type
    df["transcript_type"] = df["transcript_type"].mask(
        is_transcript_level, args.transcript_type
    )
else:
    df["gene_type"] = df["gene_type"].fillna(args.gene_type)
    df["transcript_type"] = df["transcript_type"].fillna(
        pd.Series(args.transcript_type, index=df.index).where(is_transcript_level)
    )

# Created gene and transcript features inherit the attributes of their children
df = create_missing_features(
    df,
    unique_keys={"gene": "gene_id", "transcript": "transcript_id"},
    extra_columns={"gene": GENE_ATTRIBUTES, "transcript": TRANSCRIPT_ATTRIBUTES},
    strict="raise",
)

# gtfparse reads a missing frame ('.') as 0; frame only applies to coding features
df["frame"] = df["frame"].astype(object).where(df["feature"].isin(CODING_FEATURES), ".")


# Order the features as gene, transcript, then the others within each spike-in
def sort_key(column):
    if column.name == "feature":
        return column.map(FEATURE_ORDER).fillna(len(FEATURE_ORDER))
    return column


df = df.sort_values(["seqname", "feature"], key=sort_key, kind="stable")

write_gtf(df, args.output)
