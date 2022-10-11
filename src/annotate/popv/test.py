import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'executable': './target/native/annotate/popv',
    'resources_dir': './resources_test/annotation_test_data/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna_anndata.h5ad"
ontology_files_path = f"{meta['resources_dir']}/popv_cl_ontology"

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
            "--tissue_tabula_sapiens", "Blood",
            "--tissue_reference_file", "./resources_test/annotation_test_data/TS_Blood-subsampled.h5ad",
            "--methods", "scvi,scanvi",
            "--query_obs_covariate", "batch",
            "--query_obs_cell_type_key": "none",
            "--query_obs_cell_type_unknown_label": "unknown",
            "--reference_obs_cell_type_key": "cell_ontology_class",
            "--reference_obs_covariate": "donor,method",
            "--ontology_files_path": ontology_files_path
            ])
        #self.assertTrue(Path("output.h5mu").is_file())
        #output_data = mudata.read_h5mu("output.h5mu")

if __name__ == '__main__':
    unittest.main()
