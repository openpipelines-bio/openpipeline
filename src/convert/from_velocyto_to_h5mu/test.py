import unittest
import subprocess
from pathlib import Path
import mudata
import loompy
## VIASH START
meta = {
    'functionality_name': './target/native/convert/from_velocyto_to_h5mu/from_velocyto_to_h5mu',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_loom = f"{resources_dir}/cellranger_tiny.loom"

class TestFromVelocytoToMuData(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([f"./{functionality_name}"] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_convert(self):
        self._run_and_check_output([
            "--input", input_loom,
            "--output", "foo.h5mu"])
        self.assertTrue(Path("./foo.h5mu").is_file())
        output_data = mudata.read_h5mu("./foo.h5mu")
        self.assertListEqual(list(output_data.mod.keys()), ["rna_velocity"])
        with loompy.connect(input_loom) as ds:
            self.assertTupleEqual(output_data.shape[::-1], ds.shape)

if __name__ == '__main__':
    unittest.main()