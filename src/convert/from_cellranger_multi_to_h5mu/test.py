from unittest import main, TestCase
import subprocess
from pathlib import Path
from mudata import read_h5mu

## VIASH START
meta = {
    'executable': './target/docker/convert/from_cellranger_multi_to_h5mu/from_cellranger_multi_to_h5mu',
    'resources_dir': './resources_test/10x_5k_anticmv/',
}
## VIASH END


resources_dir, executable = meta["resources_dir"], meta["executable"]
cellranger_multi_output = f"{resources_dir}/10x_5k_anticmv/processed/10x_5k_anticmv.cellranger_multi.output.output"

class TestCellrangerMultiToH5mu(TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_cellranger_multi_basic(self):
        self._run_and_check_output(["--input", cellranger_multi_output,
                                    "--output", "output.h5mu"])
        self.assertTrue(Path("output.h5mu").is_file())
        converted_data = read_h5mu("output.h5mu")
        self.assertListEqual(list(converted_data.mod.keys()), ['rna', 'prot', 'vdj_t'])
        self.assertListEqual(list(converted_data.uns.keys()), ['metrics_cellranger'])
        self.assertListEqual(converted_data.uns['metrics_cellranger'].columns.to_list(),
                             ['Category', 'Library Type', 'Grouped By', 'Group Name',
                              'Metric Name', 'Metric Value'])


        
if __name__ == "__main__":
    main()