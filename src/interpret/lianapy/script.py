import sys
import liana
import mudata

# TODO: Remove when grouping labels exist
# For sign/PCA/
import pandas as pd

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "output_compression": "gzip",
    "modality": "rna",
    "layer": None,
    "gene_symbol": "gene_symbol",
    "groupby": "harmony_integration_leiden_1.0",
    "resource_name": "consensus",
    "expr_prop": 0.1,
    "min_cells": 5,
    "aggregate_method": "rra",
    "return_all_lrs": False,
    "n_perms": 100,
}
meta = {"cpus": 4}
### VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import write_h5ad_to_h5mu_with_compression


def main():
    # Get input data
    mod = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])

    # Add dummy grouping labels when they do not exist
    if par["groupby"] not in mod.obs:
        raise ValueError(
            f"Column {par['groupy']} does not exist in "
            f".obs for modality {par['modality']}."
        )
    mod_col = mod.obs[par["groupby"]]
    original_groupby_col = mod_col.copy()
    if not isinstance(mod_col, pd.CategoricalDtype):
        mod.obs[par["groupby"]] = mod_col.astype(str).astype("category")

    # Solve gene labels
    orig_gene_label = mod.var.index
    mod.var_names = mod.var[par["gene_symbol"]].astype(str)
    mod.var_names_make_unique()

    if not meta.get("cpus"):
        meta["cpus"] = 1
    n_jobs = meta["cpus"] - 1 if meta["cpus"] > 2 else 1
    liana.mt.rank_aggregate(
        adata=mod,
        groupby=par["groupby"],
        resource_name=par["resource_name"],
        expr_prop=par["expr_prop"],
        min_cells=par["min_cells"],
        aggregate_method=par["aggregate_method"],
        consensus_opts=par["consensus_opts"],
        return_all_lrs=par["return_all_lrs"],
        layer=par["layer"],
        de_method=par["de_method"],
        n_perms=par["n_perms"],
        n_jobs=n_jobs,
        verbose=True,
        inplace=True,
        use_raw=False,
    )

    # Return original gene labels
    mod.var_names = orig_gene_label

    # Undo modifications to groupby column
    mod.obs[par["groupby"]] = original_groupby_col

    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], mod, par["output_compression"]
    )


if __name__ == "__main__":
    main()
