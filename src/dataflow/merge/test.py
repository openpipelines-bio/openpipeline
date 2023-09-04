import sys

from tempfile import NamedTemporaryFile
from pathlib import Path
import pytest
import subprocess
from pathlib import Path
import mudata as md
import pandas as pd
import numpy as np
import re

## VIASH START
meta = {
    'executable': './target/docker/dataflow/merge/merge',
    'resources_dir': './resources_test/merge_test_data/',
    'config': '/home/di/code/openpipeline/src/dataflow/merge/config.vsh.yml'
}
## VIASH END


resources_dir = meta["resources_dir"]
input_sample1_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
input_sample2_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu"

@pytest.fixture
def mudata_non_overlapping_observations(request, tmp_path):
    mudata_to_change_path = request.param
    temp_h5mu = tmp_path / "duplicate_observations.h5mu"
    mudata_to_change = md.read(mudata_to_change_path)
    # Remove 1 observation
    removed_observation_name = mudata_to_change.obs.index[-1]
    mudata_to_change = mudata_to_change[:mudata_to_change.n_obs-1]
    mudata_to_change.write(temp_h5mu, compression="gzip")
    return temp_h5mu, removed_observation_name

@pytest.fixture
def extra_var_column_value():
    return "bar"

@pytest.fixture
def mudata_with_extra_var_column(tmp_path, request, extra_var_column_value):
    [sample1_path, sample2_path] = request.param
    result = []
    for sample_path, column_value in ((sample1_path, extra_var_column_value), (sample2_path, np.nan)):
        sample = md.read_h5mu(sample_path)
        mod_names = list(sample.mod.keys())
        assert len(mod_names) == 1
        mod_name = mod_names[0]
        sample.mod[mod_name].var['test'] = column_value
        sample.var = sample.var.convert_dtypes(infer_objects=True,
                                               convert_integer=True,
                                               convert_string=False,
                                               convert_boolean=True,
                                               convert_floating=False)
        new_path = tmp_path / (Path(sample_path).stem + "_extra_col.h5mu")
        sample.write(new_path)
        result.append(new_path)
    return result


def test_merge(run_component):
    """
    Test a simple merge with fully overlapping observations
    """
    run_component([
        "--input", input_sample1_file,
        "--input", input_sample2_file,
        "--output", "merge.h5mu",
        "--output_compression", "gzip"])

    assert Path("merge.h5mu").is_file()
    concatenated_data = md.read("merge.h5mu")
    data_sample1 = md.read(input_sample1_file)
    data_sample2 = md.read(input_sample2_file)

    assert concatenated_data.n_mod == 2
    assert concatenated_data.obsm_keys() == ['prot', 'rna']
    # In this case, the observations overlap perfectly in the two input files
    assert concatenated_data.n_obs == data_sample1.n_obs
    assert concatenated_data.n_obs == data_sample2.n_obs
    pd.testing.assert_index_equal(concatenated_data.obs.index, data_sample1.obs.index)
    pd.testing.assert_index_equal(concatenated_data.obs.index, data_sample2.obs.index)

    assert set(data_sample1.var_keys()) | set(data_sample2.var_keys()) == \
           set(concatenated_data.var_keys())

    assert set(concatenated_data.var_names) ==  (set(data_sample1.var_names) | set(data_sample2.var_names))
    assert concatenated_data.var_keys() == ['gene_id', 'feature_type', 'genome']

@pytest.mark.parametrize("mudata_non_overlapping_observations", [input_sample1_file], indirect=["mudata_non_overlapping_observations"])
def test_merge_non_overlapping_observations(run_component, mudata_non_overlapping_observations):
    """
    Merge with differing observations in the samples
    """
    edited_h5mu, removed_observation_name = mudata_non_overlapping_observations
    # Remove 1 observation
    run_component([
        "--input", edited_h5mu,
        "--input", input_sample2_file,
        "--output", "merge.h5mu"])
    assert  Path("merge.h5mu").is_file()
    concatenated_data = md.read("merge.h5mu", backed=False)
    data_sample1 = md.read(edited_h5mu, backed=False)
    data_sample2 = md.read(input_sample2_file, backed=False)
    assert concatenated_data[removed_observation_name:]['rna'].n_obs == 0
    assert concatenated_data[removed_observation_name:]['prot'].n_obs == 1

    assert set(concatenated_data.obs_names) == (set(data_sample1.obs_names) | set(data_sample2.obs_names))
    np.testing.assert_equal(concatenated_data[removed_observation_name:]['rna'].X.data,
                            np.array([]))
    np.testing.assert_equal(concatenated_data.copy()[removed_observation_name:]['prot'].X.data,
                            data_sample2.copy()[removed_observation_name:]['prot'].X.data)
    
@pytest.mark.parametrize("extra_var_column_value,expected", [("bar", "bar"), (True, True), (0.1, 0.1), (np.nan, pd.NA)])
@pytest.mark.parametrize("mudata_with_extra_var_column",
                          [(input_sample1_file, input_sample2_file)],
                          indirect=["mudata_with_extra_var_column"])
def test_boolean_and_na_types(run_component, mudata_with_extra_var_column, expected):
    """
    Test if merging booleans of NAs results in the .var .obs column being writeable
    """
    input_file_1, input_file_2 = mudata_with_extra_var_column
    run_component([
        "--input", input_file_1,
        "--input", input_file_2,
        "--output", "merge.h5mu"])
    assert Path("merge.h5mu").is_file()
    merged_data = md.read("merge.h5mu", backed=False)
    first_sample_mod = list(md.read(input_file_1).mod)[0]
    second_sample_mod = list(md.read(input_file_2).mod)[0]
    if not pd.isna(expected):
        assert merged_data.var.loc['MIR1302-2HG']['test'] == expected
        assert merged_data.mod[first_sample_mod].var.loc['MIR1302-2HG']['test'] == expected
    else:
        assert pd.isna(merged_data.var.loc['MIR1302-2HG']['test'])
        assert pd.isna(merged_data.mod[first_sample_mod].var.loc['MIR1302-2HG']['test'])
    assert pd.isna(merged_data.var.loc['CD3_TotalSeqB']['test'])
    assert pd.isna(merged_data.mod[second_sample_mod].var.loc['CD3_TotalSeqB']['test'])


def test_same_modalities_raises(run_component):
    """
    Raise when trying to merge modalities with the same name.
    """
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile:
        data_sample1 = md.read(input_sample1_file)
        data_sample1.mod = {'prot': data_sample1.mod['rna']}
        data_sample1.write(tempfile.name, compression="gzip")
        with pytest.raises(subprocess.CalledProcessError) as err:
            run_component([
                "--input", tempfile.name,
                "--input", input_sample2_file,
                "--output", "merge.h5mu"])
        assert re.search(r"ValueError: Modality 'prot' was found in more than 1 sample\.",
                  err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
