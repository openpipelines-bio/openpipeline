import sys
import os
import pytest
import unittest

## VIASH START
meta = {
    'resources_dir': './resources_test/'
}
## VIASH END

OUTPUT_FILE = "output.h5mu"


def test_cellxgene(run_component):
    run_component([
        "--input_database", "CellxGene",
        "--modality", "rna",
        "--cell_ontology_release", "2023-05-22",
        "--cellxgene_release", "2023-05-15",
        "--species", "homo_sapiens",
        "--tissue", "lung",
        "--obs_column_names", "disease",
        "--output", OUTPUT_FILE,
        "--metadata_only", "True"
    ])

    # check whether file exists
    assert os.path.exists(OUTPUT_FILE), "Output file does not exist"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "--capture=no"], plugins=["viashpy"]))
