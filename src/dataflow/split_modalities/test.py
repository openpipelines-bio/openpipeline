import unittest
import subprocess
from pathlib import Path
import mudata as md
import re

## VIASH START
meta = {
    'functionality_name': './target/native/dataflow/split_modalities/split_modalities',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]

input_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


class TestSplit(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_split(self):
        output_dir = Path("foo")
        output_types = Path("foo.csv")
        self._run_and_check_output([
            "--input", input_file,
            "--output", str(output_dir),
            "--output_types", str(output_types)
        ])
        self.assertTrue(output_types.is_file())

        # todo: check whether contents of output_types is correct

        self.assertTrue(output_dir.is_dir())
        dir_content = [h5mu_file for h5mu_file in output_dir.iterdir() if h5mu_file.suffix == ".h5mu"]
        rna_file = Path("foo/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu")
        prot_file = Path("foo/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu")
        self.assertSetEqual(set(dir_content), set([prot_file, rna_file]))
        input_file_contents = md.read_h5mu(input_file)
        rna = md.read_h5mu(rna_file)
        prot = md.read_h5mu(prot_file)

        self.assertEqual(rna.n_mod, 1)
        self.assertEqual(prot.n_mod, 1)

        self.assertEqual(rna.n_obs, input_file_contents.n_obs)
        self.assertEqual(prot.n_obs, input_file_contents.n_obs)

        # When a var_key is only present for one modality, it is prefixed by the name of the
        # modality followed by a colon and the name of the key (in the global .var).
        replace_regex = r"(^rna:|^prot:)"
        self.assertSetEqual(set(rna.var_keys()) | set(prot.var_keys()),
                            set(re.sub(replace_regex, "", col_name)
                                       for col_name in input_file_contents.var_keys()))

        self.assertSetEqual(set(rna.var_keys()),
                            set(input_file_contents.mod['rna'].var.columns))
        self.assertSetEqual(set(rna.var_keys()),
                            set(input_file_contents.mod['prot'].var.columns))


if __name__ == '__main__':
    unittest.main()
