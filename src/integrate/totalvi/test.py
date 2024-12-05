import sys
import pytest
import mudata

## VIASH START
meta = {
    "executable": "./target/docker/integrate/totalvi/totalvi",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


def test_totalvi(run_component, tmp_path):
    """Map data containing proteins on itself"""
    output_path = tmp_path / "output.h5mu"
    ref_model_path = tmp_path / "totalvi_reference_model"
    query_model_path = tmp_path / "totalvi_query_model"

    run_component(
        [
            "--input",
            input_file,
            "--reference",
            input_file,
            "--query_proteins_modality",
            "prot",
            "--reference_proteins_modality",
            "prot",
            "--var_input",
            "filter_with_hvg",
            "--reference_model_path",
            str(ref_model_path),
            "--query_model_path",
            str(query_model_path),
            "--max_epochs",
            "1",
            "--max_query_epochs",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()
    output_data = mudata.read_h5mu(output_path)
    assert "X_integrated_totalvi" in output_data.mod["rna"].obsm
    assert "X_totalvi_normalized_rna" in output_data.mod["rna"].obsm
    assert "X_totalvi_normalized_protein" in output_data.mod["prot"].obsm


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
