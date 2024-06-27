import sys
import os
import pytest
import mudata as mu

## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

def test_simple_execution(run_component, random_h5mu_path):
    output_file = "output.h5mu"

    run_component([
        "--input", input_file,
        "--reference", reference_file,
        "--output", random_h5mu_path()
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
        