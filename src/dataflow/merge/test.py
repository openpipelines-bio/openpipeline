from tempfile import NamedTemporaryFile
from pathlib import Path
import unittest
import subprocess
from pathlib import Path
import mudata as md
import pandas as pd
import numpy as np

## VIASH START
meta = {
    'functionality_name': './target/docker/dataflow/merge/merge',
    'resources_dir': './resources_test/merge_test_data/'
}
## VIASH END


resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_sample1_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
input_sample2_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu"


class TestMerge(unittest.TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
                print(e.stdout.decode("utf-8"))
            raise e

    def test_merge(self):
        """
        Test a simple merge with fully overlapping observations
        """
        self._run_and_check_output([
            "--input", input_sample1_file,
            "--input", input_sample2_file,
            "--output", "merge.h5mu",
            "--output_compression", "gzip"])

        self.assertTrue(Path("merge.h5mu").is_file())
        concatenated_data = md.read("merge.h5mu")
        data_sample1 = md.read(input_sample1_file)
        data_sample2 = md.read(input_sample2_file)

        self.assertEqual(concatenated_data.n_mod, 2)
        self.assertListEqual(concatenated_data.obsm_keys(), ['prot', 'rna'])
        # In this case, the observations overlap perfectly in the two input files
        self.assertEqual(concatenated_data.n_obs, data_sample1.n_obs)
        self.assertEqual(concatenated_data.n_obs, data_sample2.n_obs)
        pd.testing.assert_index_equal(concatenated_data.obs.index, data_sample1.obs.index)
        pd.testing.assert_index_equal(concatenated_data.obs.index, data_sample2.obs.index)

        self.assertSetEqual(set(data_sample1.var_keys()) | set(data_sample2.var_keys()),
                            set(concatenated_data.var_keys()))

        self.assertSetEqual(set(concatenated_data.var_names), set(data_sample1.var_names) | set(data_sample2.var_names))
        self.assertListEqual(concatenated_data.var_keys(), ['gene_id', 'feature_type', 'genome'])

    def test_merge_non_overlapping_observations(self):
        """
        Merge with differing observations in the samples
        """
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile:
            data_sample1 = md.read(input_sample1_file)
            # Remove 1 observation
            removed_observation_name = data_sample1.obs.index[-1]
            data_sample1 = data_sample1[:data_sample1.n_obs-1]
            data_sample1.write(tempfile.name, compression="gzip")
            self._run_and_check_output([
                "--input", tempfile.name,
                "--input", input_sample2_file,
                "--output", "merge.h5mu"])
            self.assertTrue(Path("merge.h5mu").is_file())
            concatenated_data = md.read("merge.h5mu", backed=False)
            data_sample1 = md.read(tempfile.name, backed=False)
            data_sample2 = md.read(input_sample2_file, backed=False)
            self.assertEqual(concatenated_data[removed_observation_name:]['rna'].n_obs, 0)
            self.assertEqual(concatenated_data[removed_observation_name:]['prot'].n_obs, 1)
            concatenated_data[removed_observation_name:]['rna'].X.data

            self.assertSetEqual(set(concatenated_data.obs_names),
                                set(data_sample1.obs_names) | set(data_sample2.obs_names))
            np.testing.assert_equal(concatenated_data[removed_observation_name:]['rna'].X.data,
                                    np.array([]))
            np.testing.assert_equal(concatenated_data.copy()[removed_observation_name:]['prot'].X.data,
                                    data_sample2.copy()[removed_observation_name:]['prot'].X.data)

    def test_same_modalities_raises(self):
        """
        Raise when trying to merge modalities with the same name.
        """
        with NamedTemporaryFile('w', suffix=".h5mu") as tempfile:
            data_sample1 = md.read(input_sample1_file)
            data_sample1.mod = {'prot': data_sample1.mod['rna']}
            data_sample1.write(tempfile.name, compression="gzip")
            with self.assertRaises(subprocess.SubprocessError) as e:
                self._run_and_check_output([
                    "--input", tempfile.name,
                    "--input", input_sample2_file,
                    "--output", "merge.h5mu"], expected_raise=True)
            self.assertIn("ValueError: Modality 'prot' was found in more than 1 sample.",
                          e.exception.output.decode('utf-8'))

if __name__ == '__main__':
    unittest.main()
