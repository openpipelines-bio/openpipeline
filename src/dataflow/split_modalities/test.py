import sys
import pytest
import mudata as md
import re

## VIASH START
meta = {
    'functionality_name': './target/native/dataflow/split_modalities/split_modalities',
    'resources_dir': './resources_test/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_split(run_component, tmp_path):
    output_dir = tmp_path / "foo"
    output_types = tmp_path / "foo.csv"
    
    run_component([
        "--input", input_file,
        "--output", str(output_dir),
        "--output_types", str(output_types),
        "--output_compression", "gzip"
    ])
    assert output_types.is_file()
    assert output_dir.is_dir()

    # todo: check whether contents of output_types is correct

    # check output dir
    dir_content = [h5mu_file for h5mu_file in output_dir.iterdir() if h5mu_file.suffix == ".h5mu"]
    rna_file = output_dir / "pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
    prot_file = output_dir / "pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu"
    assert set(dir_content), set([prot_file == rna_file])
    input_file_contents = md.read_h5mu(input_file)
    rna = md.read_h5mu(rna_file)
    prot = md.read_h5mu(prot_file)

    assert rna.n_mod == 1
    assert prot.n_mod == 1

    assert rna.n_obs == input_file_contents.n_obs
    assert prot.n_obs == input_file_contents.n_obs

    # When a var_key is only present for one modality, it is prefixed by the name of the
    # modality followed by a colon and the name of the key (in the global .var).
    replace_regex = r"(^rna:|^prot:)"
    expected_var_keys = {re.sub(replace_regex, "", col_name) for col_name in input_file_contents.var_keys()}
    assert set(rna.var_keys()) | set(prot.var_keys()) == expected_var_keys

    assert set(rna.var_keys()) == set(input_file_contents.mod['rna'].var.columns)
    assert set(rna.var_keys()) == set(input_file_contents.mod['prot'].var.columns)

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
