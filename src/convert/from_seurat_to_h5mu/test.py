import subprocess
import sys
import pytest
import mudata as mu
import numpy as np
import re

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/executable/convert/from_h5ad_to_h5mu/from_h5ad_to_h5mu",
    "config": "./src/convert/from_h5ad_to_h5mu/config.vsh.yaml",
}
## VIASH END

ori_mdata_path = (
    meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
)
ori_seurat_path = (
    meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.rds"
)


def test_run(run_component, random_h5mu_path):
    conv_mdata_path = random_h5mu_path()

    cmd_pars = [
        "--input",
        ori_seurat_path,
        "--output",
        conv_mdata_path,
    ]
    run_component(cmd_pars)

    assert conv_mdata_path.is_file(), "No output was created."

    ori_mdata = mu.read_h5mu(ori_mdata_path)
    conv_mdata = mu.read_h5mu(conv_mdata_path)

    assert list(conv_mdata.mod.keys()) == ["rna"]

    ori_adata = ori_mdata.mod["rna"]
    conv_adata = conv_mdata.mod["rna"]

    conv_adata.obs = conv_adata.obs.drop(
        ["orig.ident", "nCount_RNA", "nFeature_RNA"], axis=1
    )

    assert all(key in list(conv_adata.layers) for key in ["X", "log_normalized"]), (
        "layers keys differ"
    )
    assert np.allclose(conv_adata.layers["X"].todense(), ori_adata.X.todense()), (
        "X layer differs"
    )
    assert conv_adata.X is None, "X layer was created"
    assert np.allclose(
        conv_adata.layers["log_normalized"].todense(),
        ori_adata.layers["log_normalized"].todense(),
    ), "log_normalized layer differs"

    assert all(ori_adata.obs.keys() == conv_adata.obs.keys()), "obs keys differ"
    obs_float_keys = [
        "scrublet_doublet_score",
        "total_counts",
        "pct_of_counts_in_top_50_vars",
    ]
    assert all(conv_adata.obs[key].dtype == float for key in obs_float_keys), (
        "obs float dtypes differ"
    )
    obs_int_keys = ["num_nonzero_vars"]
    assert all(conv_adata.obs[key].dtype == "i" for key in obs_int_keys), (
        "obs int dtypes differ"
    )
    obs_bool_keys = ["filter_with_counts", "filter_with_scrublet"]
    assert all(conv_adata.obs[key].dtype == bool for key in obs_bool_keys), (
        "obs bool dtypes differ"
    )
    obs_cat_keys = ["sample_id", "harmony_integration_leiden_1.0"]
    assert all(conv_adata.obs[key].dtype.name == "category" for key in obs_cat_keys), (
        "obs category dtypes differ"
    )

    assert all(ori_adata.var.keys() == conv_adata.var.keys()), "var keys differ"
    var_float_keys = ["obs_mean", "total_counts", "pct_dropout"]
    assert all(conv_adata.var[key].dtype == float for key in var_float_keys), (
        "var float dtypes differ"
    )
    var_int_keys = ["num_nonzero_obs"]
    assert all(conv_adata.var[key].dtype == "i" for key in var_int_keys), (
        "var int dtypes differ"
    )
    var_bool_keys = ["filter_with_counts", "filter_with_hvg"]
    assert all(conv_adata.var[key].dtype == bool for key in var_bool_keys), (
        "var bool dtypes differ"
    )
    var_cat_keys = ["gene_symbol", "feature_types", "genome"]
    assert all(conv_adata.var[key].dtype.name == "category" for key in var_cat_keys), (
        "var category dtypes differ"
    )

    assert ori_adata.obsm.keys() == conv_adata.obsm.keys(), "obsm keys differ"
    obsm_float_keys = ["X_pca", "X_umap", "X_pca_integrated", "knn_distances"]
    assert all(conv_adata.obsm[key].dtype == float for key in obsm_float_keys), (
        "obsm float dtypes differ"
    )
    obsm_int_keys = ["knn_indices"]
    assert all(conv_adata.obsm[key].dtype == "i" for key in obsm_int_keys), (
        "obsm int dtypes differ"
    )
    assert ori_adata.obsm["X_pca"].shape == conv_adata.obsm["X_pca"].shape, (
        "obsm X_pca shape differs"
    )
    assert ori_adata.obsm["X_umap"].shape == conv_adata.obsm["X_umap"].shape, (
        "obsm X_umap shape differs"
    )

    assert ori_adata.shape == conv_adata.shape, "shape differs"
    assert ori_adata.uns.keys() == conv_adata.uns.keys(), "uns keys differ"


def test_x_mapping(run_component, random_h5mu_path):
    conv_mdata_path = random_h5mu_path()

    cmd_pars = [
        "--input",
        ori_seurat_path,
        "--x_mapping",
        "X",
        "--output",
        conv_mdata_path,
    ]
    run_component(cmd_pars)

    assert conv_mdata_path.is_file(), "No output was created."

    conv_adata = mu.read_h5ad(conv_mdata_path, mod="rna")
    assert all(key in list(conv_adata.layers) for key in ["X", "log_normalized"]), (
        "layers keys differ"
    )
    assert conv_adata.X is not None, "X is None"


def test_conversion_argumens(run_component, random_h5mu_path):
    conv_mdata_path = random_h5mu_path()

    cmd_pars = [
        "--input",
        ori_seurat_path,
        "--layers_mapping",
        "False",
        "--obs_mapping",
        "False",
        "--var_mapping",
        "False",
        "--obsm_mapping",
        "False",
        "--uns_mapping",
        "False",
        "--output",
        conv_mdata_path,
    ]
    run_component(cmd_pars)

    assert conv_mdata_path.is_file(), "No output was created."

    conv_adata = mu.read_h5ad(conv_mdata_path, mod="rna")
    assert conv_adata.obs.shape[1] == 0, "obs was converted"
    assert conv_adata.var.shape[1] == 0, "var was converted"
    assert len(conv_adata.obsm) == 0, "obsm was converted"
    assert len(conv_adata.uns) == 0, "uns was converted"
    assert len(conv_adata.layers) == 0, "layers were converted"


def test_custom_mapping(run_component, random_h5mu_path):
    conv_mdata_path = random_h5mu_path()

    cmd_pars = [
        "--input",
        ori_seurat_path,
        "--layers_mapping",
        '{"layer1": "X", "layer2": "log_normalized"}',
        "--obs_mapping",
        '{"obs1": "filter_with_counts", "obs2": "filter_with_scrublet"}',
        "--var_mapping",
        '{"var1": "gene_symbol", "var2": "feature_types"}',
        "--obsm_mapping",
        '{"obsm1": "X_pca", "obsm2": "X_umap"}',
        "--uns_mapping",
        '{"uns1": "log1p", "uns2": "neighbors"}',
        "--output",
        conv_mdata_path,
    ]
    run_component(cmd_pars)

    assert conv_mdata_path.is_file(), "No output was created."

    conv_adata = mu.read_h5ad(conv_mdata_path, mod="rna")

    assert all(key in list(conv_adata.layers) for key in ["layer1", "layer2"]), (
        "layers keys differ"
    )
    assert all(key in list(conv_adata.obs) for key in ["obs1", "obs2"]), (
        "obs keys differ"
    )
    assert all(key in list(conv_adata.var) for key in ["var1", "var2"]), (
        "var keys differ"
    )
    assert all(key in list(conv_adata.obsm) for key in ["obsm1", "obsm2"]), (
        "obsm keys differ"
    )
    assert all(key in list(conv_adata.uns) for key in ["uns1", "uns2"]), (
        "uns keys differ"
    )


def test_invalid_json_mapping(run_component, random_h5mu_path):
    conv_mdata_path = random_h5mu_path()

    cmd_pars = [
        "--input",
        ori_seurat_path,
        "--layers_mapping",
        "{'layer1': 'X'}",
        "--output",
        conv_mdata_path,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)

    assert re.search(
        r"Could not parse json from argument layers_mapping",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
