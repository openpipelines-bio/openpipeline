import sys
import os
import pytest
import mudata as md


## VIASH START
meta = {
    'resources_dir': './resources_test/'
}
## VIASH END

OUTPUT_FILE = "output.h5mu"


def test_cellxgene_extract_metadata(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--tissue", "lung",
        "--obs_column_names", "disease",
        "--output", OUTPUT_FILE,
        "--metadata_only", "True"
    ])

    # check whether file exists
    assert os.path.exists(OUTPUT_FILE), "Output file does not exist"

    component_data = md.read(OUTPUT_FILE)
    obs = component_data.mod['rna'].obs
    assert "soma_joinid" in obs

def test_cellxgene_extract_metadata_expression(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--tissue", "lung",
        "--obs_column_names", "disease",
        "--output", OUTPUT_FILE,
        "--metadata_only", "False"
    ])

    # check whether file exists
    assert os.path.exists(OUTPUT_FILE), "Output file does not exist"

    component_data = md.read(OUTPUT_FILE)
    var, obs = component_data.mod['rna'].var, component_data.mod['rna'].obs
    assert "soma_joinid" in var
    assert "soma_joinid" in obs
    
if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
