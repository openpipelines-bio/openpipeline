import sys
import os
import tempfile
import shutil

import anndata
import mudata
import pandas as pd

from cellphonedb.src.core.methods import (
    cpdb_statistical_analysis_method,
    cpdb_analysis_method,
    cpdb_degs_analysis_method,
)

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "output_compression": "gzip",
    "modality": "rna",
    "layer": "log_normalized",
    "groupby": "harmony_integration_leiden_1.0",
    "gene_symbol": "gene_symbol",
    "counts_data": "hgnc_symbol",
    "method": "statistical_analysis",
    "iterations": 100,
    "threshold": 0.1,
    "pvalue": 0.05,
    "degs_file": None,
    "microenvs_file": None,
    "active_tfs_file": None,
    "score_interactions": False,
    "subsampling": False,
    "subsampling_log": False,
    "subsampling_num_pc": 100,
    "subsampling_num_cells": None,
    "output_uns_means": "cellphonedb_means",
    "output_uns_deconvoluted": "cellphonedb_deconvoluted",
    "output_uns_deconvoluted_percents": "cellphonedb_deconvoluted_percents",
    "output_uns_pvalues": "cellphonedb_pvalues",
    "output_uns_significant_means": "cellphonedb_significant_means",
    "output_uns_relevant_interactions": "cellphonedb_relevant_interactions",
    "output_uns_interaction_scores": "cellphonedb_interaction_scores",
    "output_uns_cellsign_active_interactions": "cellphonedb_cellsign_active_interactions",
    "output_uns_cellsign_active_interactions_deconvoluted": "cellphonedb_cellsign_active_interactions_deconvoluted",
}
meta = {"cpus": 1, "resources_dir": "."}
### VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import write_h5ad_to_h5mu_with_compression

CPDB_FILE_PATH = "/opt/cellphonedb/cellphonedb.zip"


def sanitize_dataframe_for_h5ad(df):
    # CellphoneDB leaves NaN in object-dtype columns such as gene_a/gene_b for
    # complex partners (which have no single gene name); h5ad can't write a
    # column that mixes strings and floats.
    df = df.copy()
    object_columns = df.select_dtypes(include="object").columns
    df[object_columns] = df[object_columns].fillna("").astype(str)
    return df


def build_counts_adata(mod):
    X = mod.layers[par["layer"]] if par["layer"] else mod.X
    gene_ids = mod.var[par["gene_symbol"]].astype(str)
    counts_adata = anndata.AnnData(
        X=X,
        obs=pd.DataFrame(index=mod.obs_names),
        var=pd.DataFrame(index=gene_ids),
    )
    counts_adata.var_names_make_unique()
    return counts_adata


def call_method(counts_adata, meta_file_path, output_dir, threads):
    common_kwargs = {
        "cpdb_file_path": CPDB_FILE_PATH,
        "meta_file_path": meta_file_path,
        "counts_file_path": counts_adata,
        "counts_data": par["counts_data"],
        "output_path": output_dir,
        "microenvs_file_path": par["microenvs_file"],
        "threshold": par["threshold"],
        "score_interactions": par["score_interactions"],
        "threads": threads,
    }

    if par["method"] == "statistical_analysis":
        result = cpdb_statistical_analysis_method.call(
            **common_kwargs,
            active_tfs_file_path=par["active_tfs_file"],
            iterations=par["iterations"],
            pvalue=par["pvalue"],
            subsampling=par["subsampling"],
            subsampling_log=par["subsampling_log"],
            subsampling_num_pc=par["subsampling_num_pc"],
            subsampling_num_cells=par["subsampling_num_cells"],
        )
        key_map = {
            "means": par["output_uns_means"],
            "deconvoluted": par["output_uns_deconvoluted"],
            "deconvoluted_percents": par["output_uns_deconvoluted_percents"],
            "pvalues": par["output_uns_pvalues"],
            "significant_means": par["output_uns_significant_means"],
        }
    elif par["method"] == "analysis":
        result = cpdb_analysis_method.call(**common_kwargs)
        key_map = {
            "means_result": par["output_uns_means"],
            "deconvoluted": par["output_uns_deconvoluted"],
            "deconvoluted_percents": par["output_uns_deconvoluted_percents"],
        }
    else:  # degs_analysis
        result = cpdb_degs_analysis_method.call(
            **common_kwargs,
            degs_file_path=par["degs_file"],
            active_tfs_file_path=par["active_tfs_file"],
        )
        key_map = {
            "means": par["output_uns_means"],
            "deconvoluted": par["output_uns_deconvoluted"],
            "deconvoluted_percents": par["output_uns_deconvoluted_percents"],
            "relevant_interactions": par["output_uns_relevant_interactions"],
            "significant_means": par["output_uns_significant_means"],
        }

    if par["score_interactions"] and "interaction_scores" in result:
        key_map["interaction_scores"] = par["output_uns_interaction_scores"]

    # CellphoneDB returns the bare `pandas.DataFrame` class (not an instance) for both
    # CellSign keys whenever it finds no active interactions to report - guard against
    # writing a class object into .uns.
    if isinstance(result.get("CellSign_active_interactions"), pd.DataFrame):
        key_map["CellSign_active_interactions"] = par[
            "output_uns_cellsign_active_interactions"
        ]
    if isinstance(
        result.get("CellSign_active_interactions_deconvoluted"), pd.DataFrame
    ):
        key_map["CellSign_active_interactions_deconvoluted"] = par[
            "output_uns_cellsign_active_interactions_deconvoluted"
        ]

    return result, key_map


def main():
    if par["method"] == "degs_analysis" and not par["degs_file"]:
        raise ValueError("--degs_file is required when --method is 'degs_analysis'.")
    if par["active_tfs_file"] and par["method"] == "analysis":
        raise ValueError(
            "--active_tfs_file (CellSign) is not supported by --method 'analysis'; "
            "use 'statistical_analysis' or 'degs_analysis'."
        )
    subsampling_args_set = (
        par["subsampling"]
        or par["subsampling_log"]
        or par["subsampling_num_pc"] != 100
        or par["subsampling_num_cells"]
    )
    if subsampling_args_set and par["method"] != "statistical_analysis":
        raise ValueError(
            "Subsampling is only supported by --method 'statistical_analysis'."
        )

    mod = mudata.read_h5ad(par["input"], mod=par["modality"])

    if par["groupby"] not in mod.obs:
        raise ValueError(
            f"Column {par['groupby']} does not exist in .obs for modality {par['modality']}."
        )

    counts_adata = build_counts_adata(mod)

    meta_df = pd.DataFrame(
        {"Cell": mod.obs_names, "cell_type": mod.obs[par["groupby"]].astype(str)}
    )

    tmp_dir = tempfile.mkdtemp()
    try:
        meta_file_path = os.path.join(tmp_dir, "meta.tsv")
        meta_df.to_csv(meta_file_path, sep="\t", index=False)

        output_dir = os.path.join(tmp_dir, "cellphonedb_output")
        os.makedirs(output_dir, exist_ok=True)

        threads = meta.get("cpus") or 1
        result, key_map = call_method(counts_adata, meta_file_path, output_dir, threads)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    for source_key, uns_key in key_map.items():
        mod.uns[uns_key] = sanitize_dataframe_for_h5ad(result[source_key])

    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], mod, par["output_compression"]
    )


if __name__ == "__main__":
    main()
