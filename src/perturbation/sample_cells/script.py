import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse
from scipy.stats import median_abs_deviation

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "input_layer": None,
    "var_gene_names": None,
    "obs_group": "perturbation_group",
    "disease_label": "disease",
    "healthy_label": "healthy",
    "skip_qc": False,
    "nmads": 5,
    "nmads_mt": 3,
    "pct_mt_max": 8.0,
    "min_genes": 5000,
    "n_disease_cells": None,
    "seed": 45,
    "batch_size": 10,
    "min_disease_cells": 1,
    "min_healthy_cells": 1,
    "output": "output.h5mu",
    "obs_output_selected": "perturbation_selected",
    "obs_output_perturb": "perturbation_perturb",
    "obs_output_batch": "perturbation_batch",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()

HB_PATTERN = r"^HB[ABDEGMQZ]\d*(?!\w)"


def get_counts(adata):
    if par["input_layer"]:
        if par["input_layer"] not in adata.layers:
            raise ValueError(
                f"Layer '{par['input_layer']}' not found in .layers. "
                f"Available: {list(adata.layers)}"
            )
        counts = adata.layers[par["input_layer"]]
    else:
        counts = adata.X
    return csr_matrix(counts) if issparse(counts) else csr_matrix(np.asarray(counts))


def get_symbols(adata):
    if par["var_gene_names"]:
        if par["var_gene_names"] not in adata.var.columns:
            raise ValueError(
                f"Column '{par['var_gene_names']}' not found in .var. "
                f"Available: {list(adata.var.columns)}"
            )
        return adata.var[par["var_gene_names"]].astype(str)
    return adata.var_names.to_series().astype(str)


def qc_metrics(counts, symbols):
    qc = ad.AnnData(
        X=counts,
        var=pd.DataFrame(
            {
                "mt": symbols.str.startswith("MT-").to_numpy(),
                "ribo": symbols.str.startswith(("RPS", "RPL")).to_numpy(),
                "hb": symbols.str.contains(HB_PATTERN, regex=True).to_numpy(),
            },
            index=pd.Index(np.arange(counts.shape[1]).astype(str)),
        ),
    )
    sc.pp.calculate_qc_metrics(
        qc, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return qc.obs


def is_outlier(values, nmads):
    mad = median_abs_deviation(values)
    centre = np.median(values)
    return (values < centre - nmads * mad) | (values > centre + nmads * mad)


def qc_pass(metrics):
    # np.logical_or.reduce instead of chained `|`: ruff would wrap the chain
    # onto lines starting with `|`, which the Nextflow runner strips (the
    # script is embedded in a Groovy string and passed through stripMargin)
    outlier = np.logical_or.reduce(
        [
            is_outlier(metrics[column].to_numpy(), par["nmads"])
            for column in [
                "log1p_total_counts",
                "log1p_n_genes_by_counts",
                "pct_counts_in_top_20_genes",
            ]
        ]
    )
    pct_mt = metrics["pct_counts_mt"].to_numpy()
    mt_outlier = is_outlier(pct_mt, par["nmads_mt"]) | (pct_mt > par["pct_mt_max"])
    return ~(outlier | mt_outlier)


def main():
    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    if par["obs_group"] not in adata.obs.columns:
        raise ValueError(
            f"Column '{par['obs_group']}' not found in .obs. Label the cells with "
            f"perturbation/label_cells first. Available: {list(adata.obs.columns)}"
        )
    labels = adata.obs[par["obs_group"]].astype(str).to_numpy()

    symbols = get_symbols(adata)
    if not par["skip_qc"] and not symbols.str.startswith("MT-").any():
        # without mitochondrial genes the mitochondrial filter passes every cell
        where = (
            f".var['{par['var_gene_names']}']"
            if par["var_gene_names"]
            else "the .var index"
        )
        raise ValueError(
            f"No mitochondrial gene ('MT-' prefix) found in {where}. "
            "Point --var_gene_names to the gene symbols, or use --skip_qc."
        )
    metrics = qc_metrics(get_counts(adata), symbols)
    if par["skip_qc"]:
        passed = np.ones(adata.n_obs, dtype=bool)
        logger.info("QC filter skipped (--skip_qc)")
    else:
        passed = qc_pass(metrics)
        logger.info("QC: %i of %i cells pass", passed.sum(), adata.n_obs)

    healthy = passed & (labels == par["healthy_label"])
    disease_qc = passed & (labels == par["disease_label"])
    n_genes = metrics["n_genes_by_counts"].to_numpy()
    disease_deep = np.flatnonzero(disease_qc & (n_genes > par["min_genes"]))
    logger.info(
        "disease: %i after QC, %i with more than %i genes; healthy: %i after QC",
        disease_qc.sum(),
        disease_deep.shape[0],
        par["min_genes"],
        healthy.sum(),
    )

    n_wanted = par["n_disease_cells"]
    if n_wanted is not None and n_wanted < disease_deep.shape[0]:
        rng = np.random.default_rng(par["seed"])
        picked = rng.permutation(disease_deep.shape[0])[:n_wanted]
        disease = np.sort(disease_deep[picked])
        logger.info("Sampled %i disease cells (seed %i)", n_wanted, par["seed"])
    else:
        disease = disease_deep
        logger.info("Keeping all %i disease cells", disease.shape[0])

    if disease.shape[0] < par["min_disease_cells"]:
        raise ValueError(
            f"Only {disease.shape[0]} disease cells are selected, "
            f"need at least {par['min_disease_cells']}."
        )
    if healthy.sum() < par["min_healthy_cells"]:
        raise ValueError(
            f"Only {healthy.sum()} healthy cells are selected, "
            f"need at least {par['min_healthy_cells']}."
        )

    perturb = np.zeros(adata.n_obs, dtype=bool)
    perturb[disease] = True

    n_batches = -(-disease.shape[0] // par["batch_size"])
    width = max(4, len(str(n_batches)))
    batches = np.full(adata.n_obs, "", dtype=object)
    for position, row in enumerate(disease):
        batches[row] = f"batch_{position // par['batch_size'] + 1:0{width}d}"
    logger.info(
        "Split %i disease cells into %i batches of at most %i",
        disease.shape[0],
        n_batches,
        par["batch_size"],
    )

    adata.obs[par["obs_output_selected"]] = healthy | perturb
    adata.obs[par["obs_output_perturb"]] = perturb
    adata.obs[par["obs_output_batch"]] = pd.Categorical(batches)

    logger.info("Writing output to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    main()
