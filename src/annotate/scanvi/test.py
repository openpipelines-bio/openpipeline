import sys
import os
import pytest
import subprocess
import re
import mudata as mu
import scanpy as sc
import anndata as ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal
## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component([
        "--input", input_file,
        "--reference", reference_file,
        "--output", output_file
    ])

    assert os.path.exists(output_file), "Output file does not exist"

    # input_mudata = mu.read_h5mu(input_file_transformed)
    # output_mudata = mu.read_h5mu(output_file)

    # assert_annotation_objects_equal(input_mudata.mod["prot"],
    #                                 output_mudata.mod["prot"])

    # assert list(output_mudata.mod["rna"].obs.keys()) == ['predicted_labels',
    #                                                      'over_clustering',
    #                                                      'majority_voting',
    #                                                      'conf_score']


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))