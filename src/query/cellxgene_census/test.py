import sys
import os
import pytest
import mudata as md
import numpy as np

## VIASH START
meta = {
    "resources_dir": "./resources_test/",
    "executable": "./target/executable/query/cellxgene_census",
    "config": "/home/di/code/openpipeline/src/query/cellxgene_census/config.vsh.yaml",
}
## VIASH END


def test_cellxgene_extract_metadata_expression(run_component, tmp_path):
    output_file = tmp_path / "output.h5mu"

    run_component(
        [
            "--obs_value_filter",
            "is_primary_data == True "
            "and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] "
            "and suspension_type == 'cell'",
            "--species",
            "homo_sapiens",
            "--add_dataset_metadata",
            "--output",
            output_file,
        ]
    )

    # check whether file exists
    assert os.path.exists(output_file), "Output file does not exist"

    mdata = md.read(output_file)

    assert "rna" in mdata.mod, "Output should contain 'rna' modality."
    assert mdata.mod["rna"].n_obs > 0, "Expected at least one cell."
    assert mdata.mod["rna"].n_vars > 0, "Expected at least one gene."

    ## check obs
    obs = mdata.mod["rna"].obs

    expected_obs = [
        "dataset_id",
        "assay",
        "assay_ontology_term_id",
        "cell_type",
        "cell_type_ontology_term_id",
        "development_stage",
        "development_stage_ontology_term_id",
        "disease",
        "disease_ontology_term_id",
        "donor_id",
        "is_primary_data",
        # "organism", "organism_ontology_term_id", # ← missing??
        "self_reported_ethnicity",
        "self_reported_ethnicity_ontology_term_id",
        "sex",
        "sex_ontology_term_id",
        "suspension_type",
        "tissue",
        "tissue_ontology_term_id",
        "tissue_general",
        "tissue_general_ontology_term_id",
        "soma_joinid",
        "collection_id",
        "collection_name",
        "collection_doi",
        "dataset_title",
    ]
    for exp_obs in expected_obs:
        assert exp_obs in obs.columns, f"Expected column '{exp_obs}' not found in .obs"

    assert np.all(obs["is_primary_data"] is True)

    ## check var
    var = mdata.mod["rna"].var
    expected_var = ["feature_id", "feature_name", "soma_joinid"]
    for exp_var in expected_var:
        assert exp_var in var.columns, f"Expected column '{exp_var}' not found in .var"

    ## check layers
    layers = mdata.mod["rna"].layers
    expected_layers = ["counts"]
    for exp_layer in expected_layers:
        assert (
            exp_layer in layers.keys()
        ), f"Expected layer '{exp_layer}' not found in .layers"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
