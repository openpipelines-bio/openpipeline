import sys
import os
import pytest
import mudata as md

## VIASH START
meta = {
    'resources_dir': './resources_test/',
    'executable': './target/docker/query/cellxgene_census',
    'config': '/home/di/code/openpipeline/src/query/cellxgene_census/config.vsh.yaml'
}
## VIASH END

OUTPUT_FILE = "output.h5mu"


def test_cellxgene_extract_metadata_expression(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--cell_query", "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616', 'CL:0002617', 'CL:0002615', 'CL:0001070', 'CL:1000310', 'CL:0000449', 'CL:1000309', 'CL:0002521', 'CL:0000448'] and assay_ontology_term_id in ['EFO:0009922', 'EFO:0009899', 'EFO:0011025', 'EFO:0030004', 'EFO:0030003'] and development_stage_ontology_term_id in ['HsapDv:0000084', 'HsapDv:0000098', 'HsapDv:0000097', 'HsapDv:0000099', 'HsapDv:0000096'] and suspension_type == 'cell'",
        "--output", OUTPUT_FILE,
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
    assert "soma_joinid" in var.columns
    assert "feature_id" in var.columns
    assert component_data.mod['rna'].n_obs

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
