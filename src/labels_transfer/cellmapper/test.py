import sys

import pytest
import re
import subprocess

import mudata
import numpy as np

## VIASH START
meta = {"resources_dir": "./resources_test/"}
## VIASH END


reference_h5mu_file = (
    f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
)
input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3/"
    "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)


@pytest.fixture
def prepared_data(tmp_path):
    embedding_key = "X_pca"

    reference_mdata = mudata.read_h5mu(reference_h5mu_file)
    reference_adata = reference_mdata.mod["rna"].copy()
    query_mdata = mudata.read_h5mu(input_file)
    query_adata = query_mdata.mod["rna"]

    reference_adata.obsm[embedding_key] = np.random.normal(
        size=(reference_adata.n_obs, 30)
    )
    query_adata.obsm[embedding_key] = np.random.normal(size=(query_adata.n_obs, 30))

    reference_mdata = mudata.MuData({"rna": reference_adata})
    query_out = tmp_path / "query.h5mu"
    ref_out = tmp_path / "reference.h5mu"
    query_mdata.write_h5mu(str(query_out))
    reference_mdata.write_h5mu(str(ref_out))

    return query_out, ref_out


def test_cellmapper_defaults(run_component, prepared_data, random_h5mu_path):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    run_component(
        [
            "--input",
            str(query_path),
            "--reference",
            str(ref_path),
            "--modality",
            "rna",
            "--input_obsm_features",
            "X_pca",
            "--reference_obsm_features",
            "X_pca",
            "--reference_obs_targets",
            "cell_type",
            "--n_neighbors",
            "5",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)
    assert "cell_type_pred" in out.mod["rna"].obs
    assert out.mod["rna"].obs["cell_type_pred"].notna().all()
    assert "cell_type_probability" in out.mod["rna"].obs
    assert out.mod["rna"].obs["cell_type_probability"].notna().all()


@pytest.mark.parametrize("kernel", ["jaccard", "gauss"])
def test_cellmapper_kernels(run_component, prepared_data, random_h5mu_path, kernel):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    run_component(
        [
            "--input",
            str(query_path),
            "--reference",
            str(ref_path),
            "--modality",
            "rna",
            "--input_obsm_features",
            "X_pca",
            "--reference_obsm_features",
            "X_pca",
            "--reference_obs_targets",
            "cell_type",
            "--n_neighbors",
            "5",
            "--kernel_method",
            kernel,
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)
    assert "cell_type_pred" in out.mod["rna"].obs


def test_cellmapper_custom_output_columns(
    run_component, prepared_data, random_h5mu_path
):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    run_component(
        [
            "--input",
            str(query_path),
            "--reference",
            str(ref_path),
            "--modality",
            "rna",
            "--input_obsm_features",
            "X_pca",
            "--reference_obsm_features",
            "X_pca",
            "--reference_obs_targets",
            "cell_type",
            "--output_obs_predictions",
            "cell_type_custom",
            "--output_obs_probability",
            "cell_type_custom_probability",
            "--n_neighbors",
            "5",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)
    obs = out.mod["rna"].obs
    assert "cell_type_custom" in obs
    assert "cell_type_custom_probability" in obs


def test_missing_reference_label(run_component, prepared_data, random_h5mu_path):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(query_path),
                "--reference",
                str(ref_path),
                "--modality",
                "rna",
                "--input_obsm_features",
                "X_pca",
                "--reference_obsm_features",
                "X_pca",
                "--reference_obs_targets",
                "non_existent_label",
                "--n_neighbors",
                "5",
                "--output",
                output,
            ]
        )

    assert re.search(
        r"Reference label 'non_existent_label' not found",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_embedding_key(run_component, prepared_data, random_h5mu_path):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(query_path),
                "--reference",
                str(ref_path),
                "--modality",
                "rna",
                "--input_obsm_features",
                "non_existent_embedding",
                "--reference_obsm_features",
                "X_pca",
                "--reference_obs_targets",
                "cell_type",
                "--n_neighbors",
                "5",
                "--output",
                output,
            ]
        )

    assert re.search(
        r"Embedding 'non_existent_embedding' not found",
        err.value.stdout.decode("utf-8"),
    )


def test_embedding_dimension_mismatch(
    run_component, prepared_data, random_h5mu_path, tmp_path
):
    query_path, ref_path = prepared_data
    output = random_h5mu_path()

    reference = mudata.read_h5mu(str(ref_path)).copy()
    reference.mod["rna"].obsm["X_pca"] = reference.mod["rna"].obsm["X_pca"][:, :-1]
    bad_ref = tmp_path / "bad_reference.h5mu"
    reference.write_h5mu(str(bad_ref))

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(query_path),
                "--reference",
                str(bad_ref),
                "--modality",
                "rna",
                "--input_obsm_features",
                "X_pca",
                "--reference_obsm_features",
                "X_pca",
                "--reference_obs_targets",
                "cell_type",
                "--n_neighbors",
                "5",
                "--output",
                output,
            ]
        )

    assert re.search(
        r"Embedding dimensions do not match",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
