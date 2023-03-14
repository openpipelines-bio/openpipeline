import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path

## VIASH START
meta = {
    'executable': './target/docker/integrate/scvi/scvi',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TestscVI(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scvi(self):
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
                "--obs_batch", "batch",
                "--output", "output.h5mu",
                "--model_output", "test/",
                "--max_epochs", "1",
                "--output_compression", "gzip"])
            self.assertTrue(Path("output.h5mu").is_file())
            output_data = mudata.read_h5mu("output.h5mu")
            self.assertIn('X_scvi_integrated', output_data.mod['rna'].obsm)
            self.assertIn('_scvi_batch', output_data.mod['rna'].obs.columns.tolist())
            self.assertIn('_scvi_labels',  output_data.mod['rna'].obs.columns.tolist())
            self.assertTrue(Path("test").is_dir())
            self.assertTrue(Path("test/model.pt").is_file())

if __name__ == '__main__':
    unittest.main()
