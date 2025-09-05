import sys
import pytest
import mudata
import subprocess
import re

## VIASH START
meta = {
    "executable": "./target/executable/integrate/scarches/scarches",
    "resources_dir": "./resources_test/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
reference = f"{meta['resources_dir']}/HLCA_reference_model.zip"
scvi_model = f"{meta['resources_dir']}/scvi_model"
scanvi_model = f"{meta['resources_dir']}/scanvi_model"


@pytest.fixture
def input_with_batch(tmp_path):
    tmp_input_path = tmp_path / "input.h5mu"

    input_data = mudata.read_h5mu(input_file)
    mod = input_data.mod["rna"]
    number_of_obs = mod.n_obs
    mod.obs["batch"] = "A"
    column_index = mod.obs.columns.get_indexer(["batch"])
    mod.obs.iloc[slice(number_of_obs // 2, None), column_index] = "B"
    input_data.write(tmp_input_path)

    return tmp_input_path, input_data


def test_hlca_reference_model(run_component, input_with_batch, tmp_path):
    tmp_input_path, _ = input_with_batch
    output_path = tmp_path / "output.h5mu"
    output_model_path = tmp_path / "model_output"

    # run component
    run_component(
        [
            "--input",
            str(tmp_input_path),
            "--input_obs_batch",
            "batch",
            "--unknown_celltype_label",
            "unlabeled",
            "--reference",
            reference,
            "--reference_class",
            "SCANVI",
            "--modality",
            "rna",
            "--output",
            str(output_path),
            "--model_output",
            str(output_model_path),
            "--max_epochs",
            "1",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = mudata.read_h5mu(output_path)
    assert "X_integrated_scanvi" in output_data.mod["rna"].obsm
    assert output_data["rna"].uns["integration_method"] == "SCANVI"
    assert (output_model_path / "model.pt").is_file()


def test_scvi_model(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    output_model_path = tmp_path / "model_output"

    # run component
    run_component(
        [
            "--input",
            input_file,
            "--input_obs_batch",
            "sample_id",
            "--reference",
            scvi_model,
            "--modality",
            "rna",
            "--output",
            str(output_path),
            "--model_output",
            str(output_model_path),
            "--max_epochs",
            "1",
            "--output_compression",
            "gzip",
            "--obsm_output",
            "X_integrated_scvi",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = mudata.read_h5mu(output_path)
    assert "X_integrated_scvi" in output_data.mod["rna"].obsm
    assert "scanvi_pred" not in output_data.mod["rna"].obs
    assert "scanvi_proba" not in output_data.mod["rna"].obs
    assert output_data["rna"].uns["integration_method"] == "SCVI"
    assert (output_model_path / "model.pt").is_file()


def test_scanvi_model(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    output_model_path = tmp_path / "model_output"

    # run component
    run_component(
        [
            "--input",
            input_file,
            "--input_obs_batch",
            "sample_id",
            "--reference",
            scanvi_model,
            "--modality",
            "rna",
            "--output",
            str(output_path),
            "--model_output",
            str(output_model_path),
            "--max_epochs",
            "1",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()

    # check output
    output_data = mudata.read_h5mu(output_path)
    assert "X_integrated_scanvi" in output_data.mod["rna"].obsm
    assert output_data["rna"].uns["integration_method"] == "SCANVI"
    assert "scanvi_pred" in output_data.mod["rna"].obs
    assert "scanvi_proba" in output_data.mod["rna"].obs

    predictions = output_data.mod["rna"].obs["scanvi_pred"]
    probabilities = output_data.mod["rna"].obs["scanvi_proba"]

    assert predictions.dtype == "category", (
        "Calculated predictions should be category dtype"
    )
    assert not all(predictions.isna()), "Not all predictions should be NA"
    assert probabilities.dtype == "float32", (
        "Calculated probabilities should be float32 dtype"
    )
    assert all(0 <= value <= 1 for value in probabilities), (
        ".obs at celltypist_probability has values outside the range [0, 1]"
    )

    assert (output_model_path / "model.pt").is_file()


def test_raises_with_missing_params(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    output_model_path = tmp_path / "model_output"

    # run component
    args = [
        "--input",
        input_file,
        "--reference",
        scvi_model,
        "--modality",
        "rna",
        "--output",
        str(output_path),
        "--model_output",
        str(output_model_path),
        "--max_epochs",
        "1",
        "--output_compression",
        "gzip",
        "--obsm_output",
        "X_integrated_scvi",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: The provided SCVI model requires `--input_obs_batch` to be provided.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
