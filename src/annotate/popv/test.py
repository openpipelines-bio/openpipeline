import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'executable': './target/native/annotate/popv',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna_anndata.h5ad"

class TestHarmonyPy(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_popv_cpu(self):
        input_data = mudata.read_h5mu(input_file)
        self._run_and_check_output([
            "--input", input_file,
            "--output", "output.h5mu",
            "--compression", "gzip",
            "--tissue", "Blood",
            "--methods", "scvi,scanvi",
            "--obs_covariate", "batch",
            "--obs_cell_type_key": "none",
            "--obs_cell_type_unknown_label": "unknown",
            "--reference_cell_type_key": "cell_ontology_class",
            "--reference_covariate": "donor,method"
            ])
        #self.assertTrue(Path("output.h5mu").is_file())
        #output_data = mudata.read_h5mu("output.h5mu")

if __name__ == '__main__':
    unittest.main()
