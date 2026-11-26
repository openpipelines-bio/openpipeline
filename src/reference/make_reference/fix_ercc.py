"""Complete an ERCC spike-in annotation so it follows the GENCODE/Ensembl layout.

ERCC annotations typically only contain one exon per spike-in with gene_id and
transcript_id. Add the matching gene and transcript features, plus gene/transcript
name and type attributes, so that stricter GTF validators (e.g. BD Rhapsody v4)
accept it.

Usage: python fix_ercc.py <input.gtf> <output.gtf>
"""

import sys

import pandas as pd
from gtfparse import create_missing_features, read_gtf, write_gtf

input_gtf, output_gtf = sys.argv[1:3]

df = read_gtf(input_gtf, result_type="pandas")
df = create_missing_features(
    df,
    unique_keys={"gene": "gene_id", "transcript": "transcript_id"},
    extra_columns={"transcript": ["gene_id"]},
    strict="raise",
)

df["gene_name"] = df["gene_id"]
df["gene_type"] = "ERCC"
is_transcript_level = df["feature"] != "gene"
df.loc[is_transcript_level, "transcript_name"] = df.loc[is_transcript_level, "gene_id"]
df.loc[is_transcript_level, "transcript_type"] = "ERCC"

# gtfparse reads a missing frame ('.') as 0; frame only applies to CDS features
df["frame"] = "."

# Order the features as gene, transcript, exon within each spike-in
df["feature"] = pd.Categorical(
    df["feature"], categories=["gene", "transcript", "exon"], ordered=True
)
df = df.sort_values(["seqname", "feature"], kind="stable")
df["feature"] = df["feature"].astype(str)

write_gtf(df, output_gtf)
