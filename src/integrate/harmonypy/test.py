import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'executable': './target/native/integrate/harmonypy/harmonypy',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TestHarmonyPy(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_harmonypy(self):
        input_data = mudata.read_h5mu(input_file)
        self._run_and_check_output([
            "--input", input_file,
            "--modality", "rna",
            "--obsm_input", "X_pca",
            "--obsm_output", "X_pca_int",
            "--obs_covariates", "leiden",
            "--output", "output.h5mu"])
        self.assertTrue(Path("output.h5mu").is_file())
        output_data = mudata.read_h5mu("output.h5mu")
        np.testing.assert_array_equal(output_data.mod['rna'].X.data, input_data.mod['rna'].X.data)
        np.testing.assert_array_equal(input_data.mod['rna'].obsm['X_pca'], output_data.mod['rna'].obsm['X_pca'])
        self.assertIn('X_pca_int', output_data.mod['rna'].obsm)
        self.assertTupleEqual(output_data.mod['rna'].obsm['X_pca_int'].shape, input_data.mod['rna'].obsm['X_pca'].shape)

if __name__ == '__main__':
    unittest.main()