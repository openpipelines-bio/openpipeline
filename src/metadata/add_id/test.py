from unittest import TestCase
import subprocess
import pandas as pd
from anndata import AnnData
from tempfile import NamedTemporaryFile
from mudata import MuData, read_h5mu
from pathlib import Path

## VIASH START
meta = {
    'executable': './target/docker/metadata/add_id/add_id',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2
}
## VIASH END

resources_dir, executable = meta["resources_dir"], meta["executable"]


class TestAddId(TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([executable] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def setUp(self) -> None:
        self.temp_h5mu = self._generate_h5mu()
        return super().setUp()

    def tearDown(self) -> None:
        self.temp_h5mu.close()
        return super().tearDown()


    def _generate_h5mu(self):
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
        obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
        var = pd.DataFrame([["a"], ["b"], ["c"]],
                           index=df.columns, columns=["Feat"])
        ad1 = AnnData(df, obs=obs, var=var)
        var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
        obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
        ad2 = AnnData(df, obs=obs2, var=var2)
        test_h5mu = NamedTemporaryFile(suffix=".h5mu")
        MuData({'mod1': ad1, 'mod2': ad2}).write_h5mu(test_h5mu.name)
        return test_h5mu

    def test_add_id(self):
        self._run_and_check_output([
                    "--input", self.temp_h5mu.name,
                    "--output", "with_id.h5mu",
                    "--input_id", "test_id",
                    "--output_compression", "gzip"
                    ])
        self.assertTrue(Path("with_id.h5mu").is_file())
        input_data = read_h5mu(self.temp_h5mu.name)
        output_data = read_h5mu("with_id.h5mu")
        self.assertIn("sample_id", output_data.obs.columns.to_list())
        self.assertSetEqual({"test_id"}, set(output_data.obs["sample_id"].to_list()))
        pd.testing.assert_index_equal(output_data.obs.index, input_data.obs.index)

    def test_add_id_obs_output(self):
        self._run_and_check_output([
                    "--input", self.temp_h5mu.name,
                    "--output", "with_id.h5mu",
                    "--input_id", "test_id",
                    "--obs_output", "test_key"
                    ])
        self.assertTrue(Path("with_id.h5mu").is_file())
        input_data = read_h5mu(self.temp_h5mu.name)
        output_data = read_h5mu("with_id.h5mu")
        self.assertIn("test_key", output_data.obs.columns.to_list())
        for mod_data in output_data.mod.values():
            self.assertIn("test_key", mod_data.obs.columns.to_list())
        self.assertSetEqual({"test_id"}, set(output_data.obs["test_key"].to_list()))
        pd.testing.assert_index_equal(output_data.obs.index, input_data.obs.index)


    def test_add_id_observations_unique(self):
        self._run_and_check_output([
                    "--input", self.temp_h5mu.name,
                    "--output", "with_id.h5mu",
                    "--input_id", "test_id",
                    "--make_observation_keys_unique"
                    ])
        self.assertTrue(Path("with_id.h5mu").is_file())
        input_data = read_h5mu(self.temp_h5mu.name)
        output_data = read_h5mu("with_id.h5mu")
        self.assertIn("sample_id", output_data.obs.columns.to_list())
        self.assertSetEqual({"test_id"}, set(output_data.obs["sample_id"].to_list()))
        pd.testing.assert_index_equal(output_data.obs.index, "test_id_" + input_data.obs.index)



