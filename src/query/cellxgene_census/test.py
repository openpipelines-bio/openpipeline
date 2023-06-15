import sys
import os
import pytest
import mudata as md
import pandas as pd

## VIASH START
meta = {
    'resources_dir': './resources_test/',
    'executable': './target/docker/query/cellxgene_census',
    'config': '/home/di/code/openpipeline/src/query/cellxgene_census/config.vsh.yaml'
}
## VIASH END

OUTPUT_FILE = "output.h5mu"


def test_cellxgene_extract_metadata(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--tissue", "mesothelial fibroblast",
        "--obs_column_names", "disease",
        "--output", OUTPUT_FILE,
        "--metadata_only", "True"
    ])

    # check whether file exists
    assert os.path.exists(OUTPUT_FILE), "Output file does not exist"

    component_data = md.read(OUTPUT_FILE)
    assert 'rna' in component_data.mod, "Output should contain 'rna' modality."
    obs = component_data.mod['rna'].obs
    var = component_data.mod['rna'].var
    assert var.empty

def test_cellxgene_extract_metadata_expression(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--tissue", "mesothelial fibroblast",
        "--obs_column_names", "disease",
        "--output", OUTPUT_FILE,
        "--metadata_only", "False"

    ])

    # check whether file exists
    assert os.path.exists(OUTPUT_FILE), "Output file does not exist"

    component_data = md.read(OUTPUT_FILE)
    assert 'rna' in component_data.mod, "Output should contain 'rna' modality."
    var, obs = component_data.mod['rna'].var, component_data.mod['rna'].obs
    assert not obs.empty, ".obs should not be empty"
    assert "is_primary_data" in obs.columns
    assert "cell_type_ontology_term_id" in obs.columns
    assert "disease" in obs.columns
    assert "soma_joinid" in var.colums
    assert "feature_id" in var.columns
    assert "cell_type_ontology_term_id" in var.columns
    assert component_data.mod['rna'].n_obs

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
