import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'functionality_name': './target/native/integrate/harmonypy/harmonypy',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_file = f"{resources_dir}/pbmc_1k_protein_v3_mms.h5mu"


class TestHarmonyPy(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([f"./{functionality_name}"] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_harmonypy(self):
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input_file:
            input_data = mudata.read_h5mu(input_file)
            mod = input_data.mod['rna']
            number_of_obs = mod.n_obs
            mod.obs['batch'] = 'A'
            column_index = mod.obs.columns.get_indexer(['batch'])
            mod.obs.iloc[slice(number_of_obs//2, None), column_index] = 'B'
            input_data.write(tempfile_input_file.name)
            self._run_and_check_output([
                "--input", tempfile_input_file.name,
                "--modality", "rna",
                "--pca_key", "X_pca",
                "--output_key_suffix", "harmonypy",
                "--sample_key", "batch",
                "--output", "output.h5mu"])
            self.assertTrue(Path("output.h5mu").is_file())
            output_data = mudata.read_h5mu("output.h5mu")
            np.testing.assert_array_equal(output_data.mod['rna'].X.data, input_data.mod['rna'].X.data)
            np.testing.assert_array_equal(input_data.mod['rna'].obsm['X_pca'], output_data.mod['rna'].obsm['X_pca'])
            self.assertIn('X_pca_harmonypy', output_data.mod['rna'].obsm)
            self.assertTupleEqual(output_data.mod['rna'].obsm['X_pca_harmonypy'].shape,
                                  input_data.mod['rna'].obsm['X_pca'].shape)

if __name__ == '__main__':
    unittest.main()