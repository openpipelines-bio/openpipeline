import unittest
import subprocess
from pathlib import Path
from mudata import read_h5mu
from tempfile import NamedTemporaryFile

## VIASH START
meta = {
    'executable': './target/docker/integrate/scanorama/scanorama',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TestScanorama(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_simple_integration(self):
        input_data = read_h5mu(input_file)
        mod = input_data.mod['rna']
        number_of_obs = mod.n_obs
        mod.obs['batch'] = 'A'
        column_index = mod.obs.columns.get_indexer(['batch'])
        mod.obs.iloc[slice(number_of_obs//2, None), column_index] = 'B'
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input_file:
            input_data.write(tempfile_input_file.name)
            self._run_and_check_output([
                "--input", tempfile_input_file.name,
                "--output", "output.h5mu",
                "--obs_batch", "batch",
                "--obsm_input", "X_pca"])
            self.assertTrue(Path("output.h5mu").is_file())
            data = read_h5mu("output.h5mu")
            self.assertTrue("X_scanorama" in data.mod['rna'].obsm)

    def test_obsm_output(self):
        input_data = read_h5mu(input_file)
        mod = input_data.mod['rna']
        number_of_obs = mod.n_obs
        mod.obs['batch'] = 'A'
        column_index = mod.obs.columns.get_indexer(['batch'])
        mod.obs.iloc[slice(number_of_obs//2, None), column_index] = 'B'
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input_file:
            input_data.write(tempfile_input_file.name)
            self._run_and_check_output([
                "--input", tempfile_input_file.name,
                "--output", "output.h5mu",
                "--obsm_output", "X_test",
                "--obs_batch", "batch",
                "--obsm_input", "X_pca"])
            self.assertTrue(Path("output.h5mu").is_file())
            data = read_h5mu("output.h5mu")
            self.assertTrue("X_test" in data.mod['rna'].obsm)

if __name__ == "__main__":
    unittest.main()