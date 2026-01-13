import mudata
import sys
import pytest
from pathlib import Path

## VIASH START
meta = {
    "executable": "./target/docker/integrate/totalvi/totalvi",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


def test_simple_execution(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    output_model_path = tmp_path / "totalvi_reference_model"

    run_component(
        [
            "--input",
            input_file,
            "--rna_modality",
            "rna",
            "--prot_modality",
            "prot",
            "--var_input",
            "filter_with_hvg",
            "--output_model",
            str(output_model_path),
            "--max_epochs",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()
    assert Path(output_model_path).is_dir()
    assert Path(f"{output_model_path}/model.pt").is_file()

    output_data = mudata.read_h5mu(output_path)
    input_data = mudata.read_h5mu(input_file)

    assert output_data.mod["rna"].n_obs == input_data.mod["rna"].n_obs, (
        "Number of observations in RNA modality does not match input data."
    )
    assert output_data.mod["rna"].n_vars == input_data.mod["rna"].n_vars, (
        "Number of variables in RNA modality does not match input data."
    )
    assert output_data.mod["prot"].n_obs == input_data.mod["prot"].n_obs, (
        "Number of observations in Protein modality does not match input data."
    )
    assert output_data.mod["prot"].n_vars == input_data.mod["prot"].n_vars, (
        "Number of variables in Protein modality does not match input data."
    )

    assert "X_integrated_totalvi" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["X_integrated_totalvi"].dtype == "f"
    assert output_data.mod["rna"].obsm["X_integrated_totalvi"].shape[1] == 20
    assert "X_totalvi_normalized_rna" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["X_totalvi_normalized_rna"].dtype == "f"
    assert "X_totalvi_normalized_protein" in output_data.mod["prot"].obsm
    assert output_data.mod["prot"].obsm["X_totalvi_normalized_protein"].dtype == "f"


def test_covariates(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--rna_modality",
            "rna",
            "--prot_modality",
            "prot",
            "--var_input",
            "filter_with_hvg",
            "--obs_categorical_covariate",
            "filter_with_scrublet",
            "--obs_continuous_covariate",
            "scrublet_doublet_score",
            "--max_epochs",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    input_data = mudata.read_h5mu(input_file)

    assert output_data.mod["rna"].n_obs == input_data.mod["rna"].n_obs, (
        "Number of observations in RNA modality does not match input data."
    )
    assert output_data.mod["rna"].n_vars == input_data.mod["rna"].n_vars, (
        "Number of variables in RNA modality does not match input data."
    )
    assert output_data.mod["prot"].n_obs == input_data.mod["prot"].n_obs, (
        "Number of observations in Protein modality does not match input data."
    )
    assert output_data.mod["prot"].n_vars == input_data.mod["prot"].n_vars, (
        "Number of variables in Protein modality does not match input data."
    )

    assert "X_integrated_totalvi" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["X_integrated_totalvi"].dtype == "f"
    assert "X_totalvi_normalized_rna" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["X_totalvi_normalized_rna"].dtype == "f"
    assert "X_totalvi_normalized_protein" in output_data.mod["prot"].obsm
    assert output_data.mod["prot"].obsm["X_totalvi_normalized_protein"].dtype == "f"


def test_parameters(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--rna_modality",
            "rna",
            "--prot_modality",
            "prot",
            "--var_input",
            "filter_with_hvg",
            "--obsm_output",
            "test_obsm_output",
            "--obsm_normalized_rna_output",
            "test_rna",
            "--obsm_normalized_protein_output",
            "test_prot",
            "--gene_dispersion",
            "gene-batch",
            "--protein_dispersion",
            "protein-batch",
            "--gene_likelihood",
            "zinb",
            "--latent_distribution",
            "ln",
            "--n_dimensions_latent_space",
            "10",
            "--max_epochs",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)

    assert "test_obsm_output" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["test_obsm_output"].dtype == "f"
    assert output_data.mod["rna"].obsm["test_obsm_output"].shape[1] == 10
    assert "test_rna" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["test_rna"].dtype == "f"
    assert "test_prot" in output_data.mod["prot"].obsm
    assert output_data.mod["prot"].obsm["test_prot"].dtype == "f"


def test_parameters(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--var_input",
            "filter_with_hvg",
            "--obsm_output",
            "test_obsm_output",
            "--obsm_normalized_rna_output",
            "test_rna",
            "--obsm_normalized_protein_output",
            "test_prot",
            "--gene_dispersion",
            "gene-batch",
            "--protein_dispersion",
            "protein-cell",
            "--gene_likelihood",
            "zinb",
            "--latent_distribution",
            "ln",
            "--n_dimensions_latent_space",
            "10",
            "--max_epochs",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)

    assert "test_obsm_output" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["test_obsm_output"].dtype == "f"
    assert output_data.mod["rna"].obsm["test_obsm_output"].shape[1] == 10
    assert "test_rna" in output_data.mod["rna"].obsm
    assert output_data.mod["rna"].obsm["test_rna"].dtype == "f"
    assert "test_prot" in output_data.mod["prot"].obsm
    assert output_data.mod["prot"].obsm["test_prot"].dtype == "f"
