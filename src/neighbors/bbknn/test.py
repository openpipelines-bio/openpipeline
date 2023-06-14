import unittest
from pathlib import Path
from mudata import read_h5mu
import subprocess
from tempfile import NamedTemporaryFile

## VIASH START
meta = {
    'executable': './target/docker/graph/bbknn/bbknn',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END
input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TestBBKNNn(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_simple_integration(self):
        input_data = read_h5mu(input_file)
        mod = input_data.mod['rna']
        if 'connectivities' in mod.obsp:
            del mod.obsp['connectivities']
        if 'distances' in mod.obsp:
            del mod.obsp['distances']
        if 'neighbors' in mod.uns:
            del mod.uns['neighbors']
        with  NamedTemporaryFile(suffix=".h5mu") as tempfile_input_file:
            input_data.write(tempfile_input_file.name)
            self._run_and_check_output([
                "--input", tempfile_input_file.name,
                "--output" ,"output.h5mu",
                "--obs_batch", "leiden",
                "--obsm_input", "X_pca",
                "--output_compression", "gzip"
            ])
            self.assertTrue(Path("output.h5mu").exists())
            data = read_h5mu("output.h5mu").mod['rna']
            self.assertTrue("connectivities" in data.obsp)
            self.assertTrue("distances" in data.obsp)
            self.assertTrue("neighbors" in data.uns)

unittest.main()