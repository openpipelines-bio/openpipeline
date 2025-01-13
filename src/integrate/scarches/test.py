import sys
import pytest
import mudata

## VIASH START
meta = {
    "executable": "./target/executable/integrate/scarches/scarches",
    "resources_dir": "./resources_test/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
reference = f"{meta['resources_dir']}/HLCA_reference_model.zip"


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
            "--reference",
            reference,
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
