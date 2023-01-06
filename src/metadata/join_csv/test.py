from unittest import main, TestCase
from tempfile import NamedTemporaryFile
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
import subprocess


## VIASH START
meta = {
    'executable': './target/docker/integrate/add_metadata/add_metadata',
}
## VIASH END

class TestAddMetadata(TestCase):
    def _run_and_check_output(self, args_as_list, expected_raise=False):
        try:
            subprocess.check_output([meta['executable']] + args_as_list,
                                    stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            if not expected_raise:
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
        obs = pd.DataFrame([["A", "sample1"], ["B", "sample2"]], index=df.index, columns=["Obs", "sample_id"])
        var = pd.DataFrame([["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
                           index=df.columns, columns=["Feat", "sample_id_var"])
        ad1 = AnnData(df, obs=obs, var=var)
        var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
        obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
        ad2 = AnnData(df, obs=obs2, var=var2)
        test_h5mu = NamedTemporaryFile(suffix=".h5mu")
        MuData({'mod1': ad1, 'mod2': ad2}).write_h5mu(test_h5mu.name)
        return test_h5mu

    def test_add_metadata_var(self):
        csv = pd.DataFrame({"id": ["sample1", "sample2"], "foo": ["v", "w"], "bar": ["x", "y"]})
        with NamedTemporaryFile(suffix=".csv") as temp_csv:
            csv.to_csv(temp_csv.name, index=False)
            self._run_and_check_output([
                    "--input", self.temp_h5mu.name,
                    "--input_csv", temp_csv.name,
                    "--output", "with_metadat.h5mu",
                    "--modality", "mod1",
                    "--var_key", "sample_id_var",
                    "--csv_key", "id"
                    ])
            result = read_h5mu("with_metadat.h5mu")
            original_data = read_h5mu(self.temp_h5mu.name)
            pd.testing.assert_frame_equal(result.mod['mod1'].var,
                                          pd.DataFrame({"Feat": ["a", "b", "c"],
                                                        "sample_id_var": ["sample1", "sample2", "sample1"],
                                                        "foo": ["v", "w", "v"],
                                                        "bar": ["x", "y", "x"]},
                                                       index=pd.Index(['var1', 'var2', 'var3']))
                                          .astype({"Feat": "object",
                                                   "sample_id_var": "category",
                                                   "foo": "category",
                                                   "bar": "category"}))
            pd.testing.assert_frame_equal(result.mod['mod1'].obs, original_data.mod['mod1'].obs)
            pd.testing.assert_frame_equal(result.mod['mod2'].obs, original_data.mod['mod2'].obs)
            pd.testing.assert_frame_equal(result.mod['mod2'].var, original_data.mod['mod2'].var)

    def test_add_metadata_matrix_sample_column(self):
        csv = pd.DataFrame({"id": ["sample1", "sample2"], "foo": ["v", "w"], "bar": ["x", "y"]})
        with NamedTemporaryFile(suffix=".csv") as temp_csv:
            csv.to_csv(temp_csv.name, index=False)
            self._run_and_check_output([
                    "--input", self.temp_h5mu.name,
                    "--input_csv", temp_csv.name,
                    "--output", "with_metadat.h5mu",
                    "--modality", "mod1",
                    "--obs_key", "sample_id",
                    "--csv_key", "id",
                    ])
            result = read_h5mu("with_metadat.h5mu")
            original_data = read_h5mu(self.temp_h5mu.name)
            pd.testing.assert_frame_equal(result.mod['mod1'].obs,
                                          pd.DataFrame({"Obs": ["A", "B"],
                                                        "sample_id": ["sample1", "sample2"],
                                                        "foo": ["v", "w"],
                                                        "bar": ["x", "y"]},
                                                       index=pd.Index(['obs1', 'obs2']))
                                          .astype({"Obs": "object",
                                                   "foo": "object",
                                                   "bar": "object"}))
            pd.testing.assert_frame_equal(result.mod['mod1'].var, original_data.mod['mod1'].var)
            pd.testing.assert_frame_equal(result.mod['mod2'].obs, original_data.mod['mod2'].obs)
            pd.testing.assert_frame_equal(result.mod['mod2'].var, original_data.mod['mod2'].var)

    def test_add_not_all_samples_in_csv_raises(self):
        csv = pd.DataFrame({"id": ["sample1", "lorem"], "foo": ["v", "w"], "bar": ["x", "y"]})
        with self.assertRaises(subprocess.CalledProcessError) as err:
            with NamedTemporaryFile(suffix=".csv") as temp_csv:
                csv.to_csv(temp_csv.name, index=False)
                out = self._run_and_check_output([
                        "--input", self.temp_h5mu.name,
                        "--input_csv", temp_csv.name,
                        "--output", "with_metadat.h5mu",
                        "--modality", "mod1",
                        "--obs_key", "sample_id",
                        "--csv_key", "id",
                        ], expected_raise=True)
        self.assertIn("Not all sample IDs selected from obs (using the column selected "
                      "with --var_key or --obs_key) were found in the csv file.",
                err.exception.stdout.decode('utf-8'))

if __name__ == "__main__":
    main()