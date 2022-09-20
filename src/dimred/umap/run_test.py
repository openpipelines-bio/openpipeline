import subprocess
import muon as mu
from unittest import TestCase, main
from pathlib import Path
from tempfile import NamedTemporaryFile
from mudata import read_h5mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

class TestPCA(TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_pca(self):
        self._run_and_check_output([
                "--input", input,
                "--output",  "output.h5mu",
                "--obsm_output", "X_foo",
                "--num_components", "26"
            ])
        self.assertTrue(Path("output.h5mu").is_file(), msg="No output was created.")
        data = mu.read_h5mu("output.h5mu")

        # check whether pca was found
        self.assertIn("X_foo", data.mod["rna"].obsm, msg="Check whether output was found in .obsm")
        self.assertTupleEqual(data.mod["rna"].obsm["X_foo"].shape, (data.n_obs, 26), msg="Check shapes")

    def test_selecting_input_layer(self):
        input_data = read_h5mu(input)
        input_data.mod['rna'].layers['test'] = input_data.mod['rna'].X
        with NamedTemporaryFile(suffix=".h5mu") as tempfile:
            input_data.write_h5mu(tempfile.name)
            self._run_and_check_output([
                    "--input", tempfile.name,
                    "--output",  "output.h5mu",
                    "--obsm_output", "test_foo",
                    "--num_components", "26",
                    "--layer", "test"
                ])
            self.assertTrue(Path("output.h5mu").is_file(), msg="No output was created.")
            data = mu.read_h5mu("output.h5mu")
            self.assertIn("test_foo", data.mod["rna"].obsm, msg="Check whether output was found in .obsm")
            self.assertTupleEqual(data.mod["rna"].obsm["test_foo"].shape, (data.n_obs, 26), msg="Check shapes")

    def test_raise_if_input_layer_is_missing(self):
        with self.assertRaises(subprocess.CalledProcessError) as err:
            self._run_and_check_output([
                    "--input", input,
                    "--output",  "output.h5mu",
                    "--obsm_output", "X_foo",
                    "--num_components", "26",
                    "--uns_neighbors", "does_not_exist"
                ], expected_raise=True)
            self.assertIn("ValueError: 'does_not_exist' was not found in .mod['rna'].uns.",
                          err.exception.stdout)

if __name__ == "__main__":
    main()