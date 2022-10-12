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

input_query_file = f"{meta['resources_dir']}/blood_test_query.h5ad"
input_reference_file = f"{meta['resources_dir']}/blood_test_reference.h5ad"
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
            "--input", input_query_file,
            "--output_dir", ".",
            "--compression", "gzip",
            "--tissue_tabula_sapiens", "",
            "--tissue_reference_file", input_reference_file,
            "--methods", "svm",
            "--query_obs_covariate", "batch",
            "--query_obs_cell_type_key": "none",
            "--query_obs_cell_type_unknown_label": "unknown",
            "--reference_obs_cell_type_key": "cell_ontology_class",
            "--reference_obs_covariate": "donor,method",
            "--ontology_files_path": ontology_files_path,
            "--plots": true
            ])
        #self.assertTrue(Path("output.h5mu").is_file())
        #output_data = mudata.read_h5mu("output.h5mu")

if __name__ == '__main__':
    unittest.main()
