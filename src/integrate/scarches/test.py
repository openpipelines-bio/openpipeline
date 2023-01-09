import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import os

## VIASH START
meta = {
    'executable': './target/docker/integrate/scarches/scarches',
    'resources_dir': './resources_test/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
reference = f"{meta['resources_dir']}/HLCA_reference_model.zip"

assert Path(input_file).is_file()
assert Path(reference).is_file()

class TestMappingToHLCA(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_hlca_reference_model(self):
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
                "--reference", reference,
                "--modality", "rna",
                "--output", "output.h5mu",
                "--model_output", "./model_output",
                "--max_epochs", "1"])
            self.assertTrue(Path("output.h5mu").is_file())
            output_data = mudata.read_h5mu("output.h5mu")
            self.assertIn('X_integrated_scanvi', output_data.mod['rna'].obsm)
            self.assertEqual(output_data["rna"].uns["integration_method"], "SCANVI")

            self.assertIn("model.pt", os.listdir("./model_output/"))

if __name__ == '__main__':
    unittest.main()
