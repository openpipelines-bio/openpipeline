import sys

import mudata as mu
import scanpy as sc

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "input_layer": None,
    "sel_clustering": "leiden_0.5",
    "output": "output.h5mu",
    "output_markers": "markers.csv",
    "output_compression": "gzip",
    "method": "t-test",
    "corr_method": "benjamini-hochberg",
    "tie_correct": False,
    "key_added": "rank_genes_groups",
    "filter_results": False,
    "min_in_group_fraction": 0.25,
    "max_out_group_fraction": 0.5,
    "min_fold_change": 1,
    "compare_abs": False,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger


logger = setup_logger()


def calculate_marker_genes(adata, par):
    input_layer = par.get("input_layer")
    if input_layer is not None:
        if input_layer not in adata.layers:
            raise ValueError(
                f"'{input_layer}' is not a layer in the input data. Available layers: {list(adata.layers.keys())}."
            )
        logger.info("Using layer '%s'", input_layer)
    else:
        logger.info("Using .X matrix")

    logger.info("Using '%s' method", par["method"])
    sc.tl.rank_genes_groups(
        adata,
        groupby=par["sel_clustering"],
        use_raw=False,
        layer=input_layer,
        method=par["method"],
        corr_method=par["corr_method"],
        tie_correct=par["tie_correct"],
        pts=True,
        key_added=par["key_added"],
    )


def filter_marker_genes(adata, par):
    logger.info("Min in-group fraction: %s", par["min_in_group_fraction"])
    logger.info("Min fold change: %s", par["min_fold_change"])
    logger.info("Max out-group fraction: %s", par["max_out_group_fraction"])
    logger.info("Compare absolute values: %s", par["compare_abs"])

    filtered_key_added = f"{par['key_added']}_filtered"
    sc.tl.filter_rank_genes_groups(
        adata,
        groupby=par["sel_clustering"],
        use_raw=False,
        key_added=filtered_key_added,
        min_in_group_fraction=par["min_in_group_fraction"],
        min_fold_change=par["min_fold_change"],
        max_out_group_fraction=par["max_out_group_fraction"],
        compare_abs=par["compare_abs"],
    )

    filtered_markers = {}
    for group in adata.obs[par["sel_clustering"]].cat.categories:
        filtered_results = sc.get.rank_genes_groups_df(
            adata, group=group, key=filtered_key_added
        )
        filtered_markers[group] = filtered_results["names"].tolist()

    filtering_params = adata.uns[filtered_key_added]["params"]
    del adata.uns[filtered_key_added]
    adata.uns[filtered_key_added] = {"params": filtering_params}
    adata.uns[filtered_key_added]["filters"] = {
        "min_in_group_fraction": par["min_in_group_fraction"],
        "min_fold_change": par["min_fold_change"],
        "max_out_group_fraction": par["max_out_group_fraction"],
        "compare_abs": par["compare_abs"],
    }
    adata.uns[filtered_key_added]["markers"] = filtered_markers

    return filtered_key_added, filtered_markers


def format_results(adata, par, filtered_key_added=None):
    results = sc.get.rank_genes_groups_df(adata, group=None, key=par["key_added"])

    if par["filter_results"]:
        results["is_marker"] = False
        for group in adata.obs[par["sel_clustering"]].cat.categories:
            group_markers = adata.uns[filtered_key_added]["markers"][group]
            is_group_marker = results["names"].isin(group_markers) & (
                results["group"] == group
            )
            results.loc[is_group_marker, "is_marker"] = True

    return results


def main(par):
    logger.info("Loading MuData from '%s'", par["input"])
    mdata = mu.read_h5mu(par["input"])

    logger.info("Selecting modality '%s'", par["modality"])
    if par["modality"] not in mdata.mod:
        raise ValueError(
            f"'{par['modality']}' is not a modality in the input data. Available modalities: {list(mdata.mod.keys())}."
        )
    adata = mdata[par["modality"]]

    logger.info("Selecting clustering column '%s'", par["sel_clustering"])
    if par["sel_clustering"] not in adata.obs.columns:
        raise ValueError(
            f"'{par['sel_clustering']}' is not a column in .obs of the input data"
        )

    logger.info("Calculating marker genes")
    calculate_marker_genes(adata, par)

    filtered_key_added = f"{par['key_added']}_filtered"
    if par["filter_results"]:
        logger.info("Filtering marker genes")
        filtered_key_added, filtered_markers = filter_marker_genes(adata, par)
        for group, markers in filtered_markers.items():
            logger.info("Group '%s': %s markers", group, len(markers))

    logger.info("Formatting results")
    results = format_results(adata, par, filtered_key_added)

    logger.info("Writing MuData output to '%s'", par["output"])
    mdata.write_h5mu(par["output"], compression=par["output_compression"])

    logger.info("Writing marker results to '%s'", par["output_markers"])
    results.to_csv(par["output_markers"], index=False)

    logger.info("Done")


if __name__ == "__main__":
    sys.exit(main(par))
