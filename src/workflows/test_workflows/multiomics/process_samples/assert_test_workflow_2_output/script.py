import sys
import mudata as mu
import pandas as pd
import pytest
from itertools import product
from operator import attrgetter
import os
import numpy as np

##VIASH START
par = {
    "input": "/home/di/code/openpipelines-multisample/work/16/def1332e1a1fd58567617645a7eb4d/merged.umap.output.h5mu",
    "orig_input": [
        "/home/di/code/openpipelines-multisample/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
        "/home/di/code/openpipelines-multisample/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    ],
}

meta = {"resources_dir": "src/base"}
sys.path.append("src/utils")
sys.path.append("src/base")
##VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def sample_1_h5mu():
    return mu.read_h5mu(par["orig_input"][0])


@pytest.fixture
def sample_2_h5mu():
    return mu.read_h5mu(par["orig_input"][1])


@pytest.fixture
def output_h5mu():
    return mu.read_h5mu(par["input"])


@pytest.fixture
def get_sample_output(request, output_h5mu, sample_1_h5mu, sample_2_h5mu):
    sample_id_to_file = {
        "pbmc": sample_1_h5mu,
        "pbmc_with_more_args": sample_2_h5mu,
    }
    modality_marker = request.node.get_closest_marker("modality")
    if modality_marker is None:
        raise ValueError("Please use @pytest.mark.modality")
    modality = modality_marker.args[0]
    layer_marker = request.node.get_closest_marker("layer")
    if layer_marker is None:
        raise ValueError("Please use @pytest.mark.layer")
    layer_name = layer_marker.args[0]
    layer_getter = (
        attrgetter("X") if layer_name == "X" else lambda ad: ad.layers[layer_name]
    )

    def _get_sample_output(sample_name):
        sample_h5mu = sample_id_to_file[sample_name]
        sample_mod = sample_h5mu[modality]
        mod_output = output_h5mu[modality]
        sample_output = mod_output[
            mod_output.obs["sample_id"] == sample_name, sample_mod.var_names
        ]
        obs_index = pd.Index(
            sample_output.obs_names.to_series().str.extract(
                f"{sample_name}_(.+)", expand=False
            )
        )
        sample_df = pd.DataFrame(
            layer_getter(sample_output).toarray(),
            index=obs_index,
            columns=sample_mod.var_names,
        )
        return sample_df

    return _get_sample_output


def test_modalities_present(output_h5mu):
    assert all(k in output_h5mu.obsm.keys() for k in ("rna", "prot"))


def test_rna_embeddings_present(output_h5mu):
    assert all(k in output_h5mu["rna"].obsm.keys() for k in ("X_pca", "X_umap"))
    assert "pca_loadings" in output_h5mu["rna"].varm.keys()
    assert "pca_variance" in output_h5mu["rna"].uns.keys()


def test_rna_layers(output_h5mu):
    assert "log_normalized" in output_h5mu["rna"].layers.keys()


@pytest.mark.layer("X")
@pytest.mark.modality("rna")
def test_rna_raw_expressions(get_sample_output, sample_1_h5mu, sample_2_h5mu):
    sample_1_df = get_sample_output("pbmc")
    sample_2_df = get_sample_output("pbmc_with_more_args")
    for input_counts, output_counts in [
        (sample_1_h5mu["rna"], sample_1_df),
        (sample_2_h5mu["rna"], sample_2_df),
    ]:
        input_counts_df = pd.DataFrame(
            input_counts.X.toarray(),
            index=input_counts.obs_names,
            columns=input_counts.var_names,
        )
        pd.testing.assert_frame_equal(input_counts_df, output_counts)


@pytest.mark.layer("log_normalized")
@pytest.mark.modality("rna")
def test_rna_log_normalized_layer(get_sample_output, sample_1_h5mu, sample_2_h5mu):
    sample_1_df = get_sample_output("pbmc")
    sample_2_df = get_sample_output("pbmc_with_more_args")
    for input_counts, output_counts in [
        (sample_1_h5mu["rna"], sample_1_df),
        (sample_2_h5mu["rna"], sample_2_df),
    ]:
        input_counts_df = pd.DataFrame(
            input_counts.X.toarray(),
            index=input_counts.obs_names,
            columns=input_counts.var_names,
        )
        input_counts_data = input_counts_df.to_numpy()
        output_counts_data = output_counts.to_numpy()
        pd.testing.assert_index_equal(input_counts_df.index, output_counts.index)
        pd.testing.assert_index_equal(input_counts_df.columns, output_counts.columns)
        # Undo log1p transformation; there are just normalized counts to compare to the original raw counts
        normalized_counts = np.expm1(output_counts_data)
        assert not np.allclose(
            input_counts_data.mean(axis=1), normalized_counts.mean(axis=1)
        )

        # When total_sum is None, after normalization, each observation (cell) has a total count equal
        # to the median of total counts for observations (cells) before normalization.
        total_sum_median = np.median(input_counts_data.sum(axis=1))
        normalized_sum = normalized_counts.sum(axis=1)
        assert np.allclose(normalized_sum, total_sum_median)

        # Check if data is log transformed
        # First normalize
        counts_per_cell = np.sum(input_counts_data, axis=1)
        counts_greater_than_zero = counts_per_cell[counts_per_cell > 0]
        median = np.median(counts_greater_than_zero)
        counts_per_cell = counts_per_cell / median
        # Broadcast array to correct dimensions
        counts_per_cell = counts_per_cell[:, None]
        # Do not allow devide by 0
        counts_per_cell = counts_per_cell.copy() + (counts_per_cell == 0)
        normalized_counts = np.true_divide(input_counts_data, counts_per_cell)
        # Check if log1p transformed normalized counts match with output
        assert np.allclose(np.log1p(normalized_counts), output_counts_data)


def test_rna_uns_slots(output_h5mu):
    assert all(
        k in output_h5mu["rna"].uns.keys()
        for k in (
            "neighbors",
            "pbmc_metrics_cellranger",
            "pbmc_with_more_args_metrics_cellranger",
            "pca_variance",
        )
    )


def test_rna_obs_slots(output_h5mu):
    obs_df = output_h5mu["rna"].obs
    mito_and_ribo_columns_prefixes = (
        "fraction_",
        "filter_",
        "total_counts_",
        "pct_",
    )

    mito_and_ribo_columns = [
        "".join(col_name_l)
        for col_name_l in product(
            mito_and_ribo_columns_prefixes, ["mitochondrial", "ribosomal"]
        )
    ]
    pct_of_counts_n = ("50", "100", "200", "500")
    pct_count_columns = [f"pct_of_counts_in_top_{i}_vars" for i in pct_of_counts_n]

    other_columns = [
        "sample_id",
        "scrublet_doublet_score",
        "filter_with_scrublet",
        "filter_with_counts",
        "num_nonzero_vars",
        "total_counts",
    ]
    assert all(
        column_name in obs_df.columns
        for column_name in pct_count_columns + mito_and_ribo_columns + other_columns
    )


def test_prot_layers_present(output_h5mu):
    assert "clr" in output_h5mu["prot"].layers.keys()


@pytest.mark.layer("X")
@pytest.mark.modality("prot")
def test_prot_raw_expressions(get_sample_output, sample_1_h5mu, sample_2_h5mu):
    sample_1_df = get_sample_output("pbmc")
    sample_2_df = get_sample_output("pbmc_with_more_args")
    for input_counts, output_counts in [
        (sample_1_h5mu["prot"], sample_1_df),
        (sample_2_h5mu["prot"], sample_2_df),
    ]:
        input_counts_df = pd.DataFrame(
            input_counts.X.toarray(),
            index=input_counts.obs_names,
            columns=input_counts.var_names,
        )
        pd.testing.assert_frame_equal(input_counts_df, output_counts)


@pytest.mark.layer("clr")
@pytest.mark.modality("prot")
def test_prot_log_normalized_expressions(
    get_sample_output, sample_1_h5mu, sample_2_h5mu
):
    sample_1_df = get_sample_output("pbmc")
    sample_2_df = get_sample_output("pbmc_with_more_args")
    for input_counts, output_counts in [
        (sample_1_h5mu["prot"], sample_1_df),
        (sample_2_h5mu["prot"], sample_2_df),
    ]:
        input_counts_df = pd.DataFrame(
            input_counts.X.toarray(),
            index=input_counts.obs_names,
            columns=input_counts.var_names,
        )
        pd.testing.assert_index_equal(input_counts_df.index, output_counts.index)
        pd.testing.assert_index_equal(input_counts_df.columns, output_counts.columns)
        input_counts_df_values = input_counts_df.to_numpy()
        clr_transformed = np.log1p(
            input_counts_df_values
            / np.exp(
                np.log1p(input_counts_df_values).sum(axis=0, keepdims=True)
                / input_counts_df_values.shape[0]
            )
        )
        np.testing.assert_allclose(
            clr_transformed, output_counts, rtol=1e-03, atol=1e-05
        )


if __name__ == "__main__":
    sys.exit(
        pytest.main(
            [
                "--import-mode=importlib",
                "-o",
                "python_files=script*.py .viash_script.py",
                os.path.dirname(__file__),
            ]
        )
    )
