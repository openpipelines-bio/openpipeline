import mudata as mu
import tempfile
import subprocess
import os
import sys
import numpy as np
from scipy.sparse import csr_matrix

## VIASH START
import muon

file_raw = (
    "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5"
)
mdat = muon.read_10x_h5(file_raw)
mdat = mdat[0:100000,]  # subsample to reduce computational time
file_input = "cellbender_remove_background_input.h5mu"
mdat.write_h5mu(file_input)

par = {
    # inputs
    "input": file_input,
    "modality": "rna",
    # outputs
    "output": "output.h5mu",
    "layer_output": "corrected",
    "obs_latent_rt_efficiency": "latent_rt_efficiency",
    "obs_latent_cell_probability": "latent_cell_probability",
    "obs_latent_scale": "latent_scale",
    "var_ambient_expression": "ambient_expression",
    # "obsm_latent_gene_encoding": "latent_gene_encoding",
    # args
    "total_droplets_included": None,
    "min_counts": 1000,
    "epochs": 5,
    "fpr": 0.01,
    "exclude_antibody_capture": False,
    "learning_rate": 0.001,
    "layer_corrected": "corrected",
    "cuda": False,
    "expected_cells": None,
    "model": "full",
    "low_count_threshold": 15,
    "z_dim": 100,
    "z_layers": [500],
    "training_fraction": 0.9,
    "empty_drop_training_fraction": 0.5,
    "expected_cells_from_qc": True,
    "output_compression": "gzip",
    "obsm_latent_gene_encoding": "cellbender_latent_gene_encoding",
}
meta = {
    "temp_dir": os.getenv("VIASH_TEMP"),
    "resources_dir": "src/correction/cellbender_remove_background",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

from helper import anndata_from_h5

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Performing log transformation on modality %s", mod)
data = mdata.mod[mod]

# with pathlib.Path(meta["temp_dir"]) / "cellbender" as temp_dir:
#     os.mkdir(temp_dir)
with tempfile.TemporaryDirectory(
    prefix="cellbender-", dir=meta["temp_dir"]
) as temp_dir:
    # construct paths within tempdir
    input_file = os.path.join(temp_dir, "input.h5ad")
    output_file = os.path.join(temp_dir, "output.h5")

    logger.info("Creating AnnData input file for CellBender: '%s'", input_file)
    data.write_h5ad(input_file)

    logger.info("Constructing CellBender command")
    cmd_pars = [
        "cellbender",
        "remove-background",
        "--input",
        input_file,
        "--output",
        output_file,
    ]

    extra_args = [
        ("--expected-cells", "expected_cells", True),
        ("--total-droplets-included", "total_droplets_included", True),
        ("--model", "model", True),
        ("--epochs", "epochs", True),
        ("--cuda", "cuda", False),
        ("--low-count-threshold", "low_count_threshold", True),
        ("--z-dim", "z_dim", True),
        ("--z-layers", "z_layers", True),
        ("--training-fraction", "training_fraction", True),
        ("--exclude-antibody-capture", "exclude_antibody_capture", False),
        ("--learning-rate", "learning_rate", True),
        ("--empty-drop-training-fraction", "empty_drop_training_fraction", True),
    ]
    for flag, name, is_kwarg in extra_args:
        if par[name]:
            values = par[name] if isinstance(par[name], list) else [par[name]]
            cmd_pars += [flag] + [str(val) for val in values] if is_kwarg else [flag]

    if par["expected_cells_from_qc"] and "metrics_cellranger" in data.uns:
        assert par["expected_cells"] is None, (
            "If min_counts is defined, expected_cells should be undefined"
        )
        assert par["total_droplets_included"] is None, (
            "If min_counts is defined, expected_cells should be undefined"
        )
        met = data.uns["metrics_cellranger"]
        col_name = "Estimated Number of Cells"
        assert col_name in met.columns, (
            "%s should be a column in .obs[metrics_cellranger]"
        )
        est_cells = met[col_name].values[0]
        logger.info(
            "Selecting --expected-cells %d and --total-droplets-included %d",
            est_cells,
            est_cells * 5,
        )
        cmd_pars += [
            "--expected-cells",
            str(est_cells),
            "--total-droplets-included",
            str(5 * est_cells),
        ]

    logger.info("Running CellBender: '%s'", " ".join(cmd_pars))
    out = subprocess.check_output(cmd_pars).decode("utf-8")

    logger.info("Reading CellBender 10xh5 output file: '%s'", output_file)
    # have to use custom read_10x_h5 function for now
    # will be fixed when https://github.com/scverse/scanpy/pull/2344 is merged
    # adata_out = sc.read_10x_h5(output_file, gex_only=False)
    adata_out = anndata_from_h5(output_file, analyzed_barcodes_only=False)

    logger.info("Copying X output to MuData")
    data.layers[par["layer_output"]] = adata_out.X

    logger.info("Copying .obs output to MuData")
    obs_store = {
        "obs_latent_rt_efficiency": "latent_RT_efficiency",
        "obs_latent_cell_probability": "latent_cell_probability",
        "obs_latent_scale": "latent_scale",
    }
    for to_name, from_name in obs_store.items():
        if par[to_name]:
            if from_name in adata_out.obs:
                data.obs[par[to_name]] = adata_out.obs[from_name]
            # when using unfiltered data, the values will be in uns instead of obs
            elif (
                from_name in adata_out.uns
                and "barcode_indices_for_latents" in adata_out.uns
            ):
                vec = np.zeros(data.n_obs)
                vec[adata_out.uns["barcode_indices_for_latents"]] = adata_out.uns[
                    from_name
                ]
                data.obs[par[to_name]] = vec

    logger.info("Copying .var output to MuData")
    var_store = {"var_ambient_expression": "ambient_expression"}
    for to_name, from_name in var_store.items():
        if par[to_name]:
            data.var[par[to_name]] = adata_out.var[from_name]

    logger.info("Copying obsm_latent_gene_encoding output to MuData")
    obsm_store = {"obsm_latent_gene_encoding": "latent_gene_encoding"}
    for to_name, from_name in obsm_store.items():
        if par[to_name]:
            if from_name in adata_out.obsm:
                data.obsm[par[to_name]] = adata_out.obsm[from_name]
            elif (
                from_name in adata_out.uns
                and "barcode_indices_for_latents" in adata_out.uns
            ):
                matrix_to_store = adata_out.uns[from_name]
                number_of_obs = data.X.shape[0]
                latent_space_sparse = csr_matrix(
                    (number_of_obs, par["z_dim"]), dtype=adata_out.uns[from_name].dtype
                )
                obs_rows_in_space_representation = adata_out.uns[
                    "barcode_indices_for_latents"
                ]
                latent_space_sparse[obs_rows_in_space_representation] = adata_out.uns[
                    from_name
                ]
                data.obsm[par[to_name]] = latent_space_sparse
            else:
                raise RuntimeError(
                    "Requested to save latent gene encoding, but the data is either missing "
                    "from cellbender output or in an incorrect format."
                )


logger.info("Writing to file %s", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
