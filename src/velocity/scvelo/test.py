import unittest
import subprocess
from pathlib import Path
from mudata import read_h5mu

## VIASH START
meta = {
    'functionality_name': './target/native/projection/scvelo/scvelo',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_loom = f"{resources_dir}/cellranger_tiny.loom"

class TestScVelo(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scvelo(self):
        self._run_and_check_output([
            "--input", input_loom,
            "--output", "./foo/",
            "--output_compression", "gzip"])
        self.assertTrue(Path("./foo").is_dir())
        self.assertTrue(Path("./foo/scvelo_proportions.pdf").is_file())
        self.assertTrue(Path("./foo/scvelo_embedding.pdf").is_file())
        self.assertTrue(Path("./foo/scvelo_graph.pdf").is_file())
        self.assertTrue(Path("./foo/proportions.txt").is_file())
        self.assertTrue(Path("./foo/foo.h5mu").is_file())
        
        output_data = read_h5mu("./foo/foo.h5mu")
        self.assertTrue("rna_velocity" in output_data.mod.keys())
        
if __name__ == '__main__':
    unittest.main()