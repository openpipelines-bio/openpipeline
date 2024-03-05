import sys
import pytest
from pathlib import Path
import subprocess
from mudata import MuData, read_h5mu
import pandas as pd
import numpy as np
import re
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import os

## VIASH START
meta = {
    'executable': './target/docker/dataflow/merge/merge',
    'resources_dir': './resources_test/merge_test_data/',
    'config': './src/dataflow/merge/config.vsh.yml'
}
## VIASH END

# TODO: Use synthetic data
input_sample1_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
input_sample2_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu"

@pytest.fixture
def mudata_non_overlapping_observations(request, random_h5mu_path):
    mudata_to_change_path = request.param
    temp_h5mu = random_h5mu_path()
    mudata_to_change = read_h5mu(mudata_to_change_path)
    # Remove 1 observation
    removed_observation_name = mudata_to_change.obs.index[-1]
    mudata_to_change = mudata_to_change[:mudata_to_change.n_obs-1]
    mudata_to_change.write(temp_h5mu, compression="gzip")
    return temp_h5mu, removed_observation_name

@pytest.fixture
def extra_var_column_value():
    return "bar"

@pytest.fixture
def mudata_with_extra_var_column(random_h5mu_path, request, extra_var_column_value):
    [sample1_path, sample2_path] = request.param
    result = []
    for sample_path, column_value in ((sample1_path, extra_var_column_value), (sample2_path, np.nan)):
        sample = read_h5mu(sample_path)
        mod_names = list(sample.mod.keys())
        assert len(mod_names) == 1
        mod_name = mod_names[0]
        sample.mod[mod_name].var['test'] = column_value
        sample.var = sample.var.convert_dtypes(infer_objects=True,
                                               convert_integer=True,
                                               convert_string=False,
                                               convert_boolean=True,
                                               convert_floating=False)
        new_path = random_h5mu_path()
        sample.write(new_path)
        result.append(new_path)
    return result




def test_merge(run_component, random_h5mu_path, split_small_mudata_path):
    """
    Test a simple merge with fully overlapping observations
    """
    
    output_path = random_h5mu_path()
    input_sample1_file, input_sample2_file = split_small_mudata_path
    args = [
        "--input", input_sample1_file,
        "--input", input_sample2_file,
        "--output", output_path,
        "--output_compression", "gzip"
    ]
    run_component(args)

    assert output_path.is_file()
    concatenated_data = read_h5mu(output_path)
    data_sample1 = read_h5mu(input_sample1_file)
    data_sample2 = read_h5mu(input_sample2_file)
    
    expected_concatenated_data = MuData({'mod1': data_sample1.mod['mod1'], 'mod2': data_sample2.mod['mod2']})
    assert_annotation_objects_equal(concatenated_data, expected_concatenated_data)

 
@pytest.mark.parametrize("mudata_non_overlapping_observations", [input_sample1_file], indirect=["mudata_non_overlapping_observations"])
def test_merge_non_overlapping_observations(run_component, mudata_non_overlapping_observations, random_h5mu_path):
    """
    Merge with differing observations in the samples
    """
    edited_h5mu, removed_observation_name = mudata_non_overlapping_observations
    output_path = random_h5mu_path()
    # Remove 1 observation
    run_component([
        "--input", edited_h5mu,
        "--input", input_sample2_file,
        "--output", output_path])
    
    assert  output_path.is_file()
    concatenated_data = read_h5mu(output_path, backed=False)
    data_sample1 = read_h5mu(edited_h5mu, backed=False)
    data_sample2 = read_h5mu(input_sample2_file, backed=False)
    
    expected_concatenated_data = MuData({'rna': data_sample1.mod['rna'], 'prot': data_sample2.mod['prot']})
    
    assert_annotation_objects_equal(concatenated_data, expected_concatenated_data)
    assert concatenated_data[removed_observation_name:]['rna'].n_obs == 0
    assert concatenated_data[removed_observation_name:]['prot'].n_obs == 1
    
    
@pytest.mark.parametrize("extra_var_column_value,expected", [("bar", "bar"), (True, True), (0.1, 0.1), (np.nan, pd.NA)])
@pytest.mark.parametrize("mudata_with_extra_var_column",
                          [(input_sample1_file, input_sample2_file)],
                          indirect=["mudata_with_extra_var_column"])
def test_boolean_and_na_types(run_component, mudata_with_extra_var_column, expected, random_h5mu_path):
    """
    Test if merging booleans of NAs results in the .var .obs column being writeable
    """
    input_file_1, input_file_2 = mudata_with_extra_var_column
    output_path = random_h5mu_path()

    run_component([
        "--input", input_file_1,
        "--input", input_file_2,
        "--output", output_path])
    assert output_path.is_file()
    merged_data = read_h5mu(output_path, backed=False)
    first_sample_mod = list(read_h5mu(input_file_1).mod)[0]
    second_sample_mod = list(read_h5mu(input_file_2).mod)[0]
    
    expected_merged_data = MuData({'rna': read_h5mu(input_file_1).mod['rna'],
                                   'prot': read_h5mu(input_file_2).mod['prot']})
    
    if not pd.isna(expected):
        assert merged_data.var.loc['MIR1302-2HG']['test'] == expected
        assert merged_data.mod[first_sample_mod].var.loc['MIR1302-2HG']['test'] == expected
    else:
        assert pd.isna(merged_data.var.loc['MIR1302-2HG']['test'])
        assert pd.isna(merged_data.mod[first_sample_mod].var.loc['MIR1302-2HG']['test'])
    assert pd.isna(merged_data.var.loc['CD3_TotalSeqB']['test'])
    assert pd.isna(merged_data.mod[second_sample_mod].var.loc['CD3_TotalSeqB']['test'])

    assert_annotation_objects_equal(merged_data, expected_merged_data)
    

def test_same_modalities_raises(run_component, random_h5mu_path):
    """
    Raise when trying to merge modalities with the same name.
    """
    input_sample1_tmp_path = random_h5mu_path()
    output_path = random_h5mu_path()
    data_sample1 = read_h5mu(input_sample1_file)
    data_sample1.mod = {'prot': data_sample1.mod['rna']}
    data_sample1.write(input_sample1_tmp_path.name, compression="gzip")
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_sample1_tmp_path.name,
            "--input", input_sample2_file,
            "--output", output_path])
    assert re.search(r"ValueError: Modality 'prot' was found in more than 1 sample\.",
                err.value.stdout.decode('utf-8'))
    
    if os.path.exists(input_sample1_tmp_path.name):
        os.remove(input_sample1_tmp_path.name)

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
