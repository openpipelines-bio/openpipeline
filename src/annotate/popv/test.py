import unittest
import subprocess
import scanpy as sc
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/annotate/popv/popv',
    'resources_dir': './resources_test/annotation_test_data/'
}
## VIASH END

input_query_file = "{}/annotation_test_data/blood_test_query.h5ad".format(meta['resources_dir'])
input_reference_file = "{}/annotation_test_data/blood_test_reference.h5ad".format(meta['resources_dir'])
ontology_files_path = "{}/popv_cl_ontology/".format(meta['resources_dir'])

class TestPopV(unittest.TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_popv_cpu(self):
        self._run_and_check_output([
            "--input", input_query_file,
            "--output_dir", "output/",
            "--compression", "gzip",
            "--tissue_tabula_sapiens", "",
            "--tissue_reference_file", input_reference_file,
            "--methods", "svm",
            "--query_obs_covariates", "batch",
            "--query_obs_cell_type_key", "none",
            "--query_obs_cell_type_unknown_label", "unknown",
            "--reference_obs_cell_type_key", "cell_ontology_class",
            "--reference_obs_covariate", "donor,method",
            "--ontology_files_path", ontology_files_path,
            "--plots", 'True'
            ])
        
        # check for cell typed data
        self.assertTrue(Path("output/blood_test_query_cell_typed.h5ad").is_file())

        # check confusion matrix
        self.assertTrue(Path("output/confusion_matrices.pdf").is_file())
        
        # check for per cell type agreement barplot
        self.assertTrue(Path("output/percelltype_agreement_barplot.pdf").is_file())

        # check for prediction score barplot
        self.assertTrue(Path("output/prediction_score_barplot.pdf").is_file())
        
        # check for predictions
        self.assertTrue(Path("output/predictions.csv").is_file())
        
        input_data = sc.read_h5ad(Path("output/blood_test_query_cell_typed.h5ad"))
        self.assertIn('popv_prediction', input_data.obs)
        
        

if __name__ == '__main__':
    unittest.main()
