from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
import unittest

import mudata
import numpy as np
import pandas as pd

## VIASH START
meta = {
    'executable': './target/docker/integrate/scvi_tools/scvi_tools',
    'resources_dir': './resources_test/scvi_tools/'
}
## VIASH END

query = f"{meta['resources_dir']}/pancreas_ref_test.h5mu"
reference = f"{meta['resources_dir']}/pancreas_query_test.h5mu"

class TestscVI(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scvi(self):    
        self._run_and_check_output([
            "--query", query,
            "--reference", reference,
            "--modality", "rna",
            "--obs_batch", "tech",
            "--output", "output.h5mu",
            "--max_epochs", "1"])
        
        self.assertTrue(Path("output.h5mu").is_file())
        output_data = mudata.read_h5mu("output.h5mu")
        self.assertIn('X_scvi_integrated', output_data.mod['rna'].obsm)
        self.assertIn('_scvi_batch', output_data.mod['rna'].obs.columns.tolist())
        self.assertIn('_scvi_labels',  output_data.mod['rna'].obs.columns.tolist())


class TestscanVI(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scanvi(self):    
        self._run_and_check_output([
            "--query", query,
            "--reference", reference,
            "--base_model", "scanvi",
            "--labels_key", "celltype",
            "--modality", "rna",
            "--obs_batch", "tech",
            "--output", "output.h5mu",
            "--max_epochs", "1"])
        
        self.assertTrue(Path("output.h5mu").is_file())
        output_data = mudata.read_h5mu("output.h5mu")
        self.assertIn('X_scanvi_integrated', output_data.mod['rna'].obsm)
        self.assertIn('_scvi_batch', output_data.mod['rna'].obs.columns.tolist())
        self.assertIn('_scvi_labels',  output_data.mod['rna'].obs.columns.tolist())
        self.assertIn('predicted_labels',  output_data.mod['rna'].obs.columns.tolist())


class TesttotalVI(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_totalvi(self):   
        # Create some fake protein data
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_ref:
            ref_data = mudata.read_h5mu(reference)
            protein_data = pd.DataFrame(
                columns=["a", "b", "c"],
                index=ref_data.mod["rna"].obs_names,
                data=np.random.randint(1, 10, size=(ref_data.n_obs, 3))
            )

            ref_data.mod["rna"].obsm["protein_expression"] = protein_data
            mudata.write(tempfile_ref.name, ref_data)

            self._run_and_check_output([
                "--query", query,
                "--reference", tempfile_ref.name,
                "--modality", "rna",
                "--obs_batch", "batch",
                "--output", "output.h5mu",
                "--query_proteins_key", "protein_expression",
                "--max_epochs", "1"])
            
            self.assertTrue(Path("output.h5mu").is_file())
            output_data = mudata.read_h5mu("output.h5mu")
            self.assertIn('X_scvi_integrated', output_data.mod['rna'].obsm)
            self.assertIn('_scvi_batch', output_data.mod['rna'].obs.columns.tolist())
            self.assertIn('_scvi_labels',  output_data.mod['rna'].obs.columns.tolist())


if __name__ == '__main__':
    unittest.main()
