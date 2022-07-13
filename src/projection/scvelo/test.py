import unittest
import subprocess

## VIASH START
meta = {
    'functionality_name': './target/native/projection/scvelo/scvelo',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_loom = f"{resources_dir}/rna_velocity/velocyto_processed/cellranger_tiny.loom"

class TestScVelo(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([f"./{functionality_name}"] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scvelo(self):
        self._run_and_check_output([
            "--input", input_loom,
            "--output", "./foo/"])
        raise NotImplementedError

if __name__ == '__main__':
    unittest.main()