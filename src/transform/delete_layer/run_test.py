from unittest import TestCase, main
from mudata import read_h5mu
from tempfile import NamedTemporaryFile
from pathlib import Path
from subprocess import check_output, CalledProcessError, STDOUT

## VIASH START
meta = {
    'functionality_name': './target/native/transform/delete_layer/delete_layer',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_file = f"{resources_dir}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

class TestDeleteLayer(TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            check_output([meta['executable']] + args_as_list, stderr=STDOUT)
        except CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_delete_layer(self):
        original_input_data = read_h5mu(input_file)
        new_layer = original_input_data.mod['rna'].X
        original_input_data.mod['rna'].layers['test'] = new_layer
        with NamedTemporaryFile(suffix=".h5mu") as tempfile:
            original_input_data.write_h5mu(tempfile.name)
            self._run_and_check_output([
                "--input", tempfile.name,
                "--modality", "rna",
                "--layer", "test",
                "--output", "deleted_layer.h5mu",
                "--output_compression", "gzip"])
        self.assertTrue(Path("deleted_layer.h5mu").is_file())
        output_data = read_h5mu('deleted_layer.h5mu')
        self.assertNotIn('test', output_data.mod['rna'].layers.keys())

    def test_missing_layer_raises(self):
        with self.assertRaises(CalledProcessError) as err:
            self._run_and_check_output([
                "--input", input_file,
                "--modality", "rna",
                "--layer", "test",
                "--output", "missing_layer.h5mu"],
                expected_raise=True)
            self.assertIn("Layer 'test' is not present in modality rna.",
                          err.exception.stdout.decode('utf-8'))

    def test_missing_layer_missing_ok(self):
        self._run_and_check_output([
            "--input", input_file,
            "--modality", "rna",
            "--layer", "test",
            "--output", "missing_layer_ok.h5mu",
            "--missing_ok"])
        self.assertTrue(Path("missing_layer_ok.h5mu").is_file())
        output_data = read_h5mu('missing_layer_ok.h5mu')
        self.assertNotIn('test', output_data.mod['rna'].layers.keys())


if __name__ == "__main__":
    main()