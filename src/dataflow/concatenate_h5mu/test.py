import mudata as md
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import pytest
import re
import sys
import uuid
import muon

## VIASH START
meta = {
    'executable': './target/docker/dataflow/concatenate_h5mu/concatenate_h5mu',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2,
    'config': './src/dataflow/concatenate_h5mu/config.vsh.yaml'
}
## VIASH END

meta['cpus'] = 1 if not meta['cpus'] else meta['cpus']

# Note: the .var for these samples have no overlap, so there are no conflicting annotations
# for the features that need to be handled by the concat component.
# The tests below that specifically test the concatenation of conflicting data need to introduce
# the conflict.
input_sample1_file = f"{meta['resources_dir']}/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"
input_sample2_file = f"{meta['resources_dir']}/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"

@pytest.fixture
def anndata_to_sparse_dataframe():
    def wrapper(anndata_object):
        return pd.DataFrame.sparse.from_spmatrix(anndata_object.X,
                                                 index=anndata_object.obs_names, 
                                                 columns=anndata_object.var_names)
    return wrapper

@pytest.fixture
def mudata_without_genome(tmp_path, request):
    mudatas_to_change, modalities = request.param
    result = []
    for mudata_to_change in mudatas_to_change:
        new_mudata = md.read(mudata_to_change)
        for mod in modalities:
            new_mudata.mod[mod].var.drop('genome', axis=1, inplace=True)
        # Note: when a column is present in the feature annotation for one
        # modality (MuData.mod[...].var) but not for another other, then the global
        # var (MuData.var) gets a new column 'mod_name:column_name'
        # (here atac:genome) next to the old 'column_name' (here just 'genome')
        new_mudata.update_var()
        new_mudata.var.drop('genome', axis=1, inplace=True)
        new_mudata = new_mudata[0:500,] # subsample to reduce memory consumption
        new_path = tmp_path / Path(mudata_to_change).name
        new_mudata.write(new_path, compression="gzip")
        result.append(new_path)
    return result


@pytest.fixture
def make_obs_names_unique():
    def wrapper(mudata):
        for mod_data in mudata.mod.values():
            mod_data.obs.index = mod_data.obs.index.map(
                lambda index_val: uuid.uuid4().hex + index_val
            )
    return wrapper

@pytest.fixture
def mudata_copy_with_unique_obs(request, make_obs_names_unique):
    mudata_to_copy = request.param
    mudata_contents = md.read(mudata_to_copy)
    copied_contents = mudata_contents.copy()
    make_obs_names_unique(copied_contents)
    return mudata_contents, copied_contents

@pytest.fixture
def extra_column_value_sample1():
    return "bar"

@pytest.fixture
def extra_column_value_sample2():
    return "bar"

@pytest.fixture
def extra_column_annotation_matrix():
    return 'var'

@pytest.fixture
def copied_mudata_with_extra_annotation_column(tmp_path, mudata_copy_with_unique_obs,
                                               extra_column_annotation_matrix,
                                               extra_column_value_sample1,
                                               extra_column_value_sample2):
    [original_mudata, copied_mudata] = mudata_copy_with_unique_obs
    getattr(original_mudata.mod['rna'], extra_column_annotation_matrix)['test']  = np.nan
    getattr(original_mudata.mod['atac'], extra_column_annotation_matrix)['test'] = extra_column_value_sample1
    original_mudata.update()

    getattr(copied_mudata.mod['rna'], extra_column_annotation_matrix)['test'] = np.nan
    getattr(copied_mudata.mod['atac'], extra_column_annotation_matrix)['test'] = extra_column_value_sample2
    mudata_cast_types = getattr(original_mudata, extra_column_annotation_matrix).convert_dtypes(
        infer_objects=True,
        convert_integer=True,
        convert_string=False,
        convert_boolean=True,
        convert_floating=False
    )
    setattr(original_mudata, extra_column_annotation_matrix, mudata_cast_types)
    copied_mudata.update()

    original_mudata_path = tmp_path / "original.h5mu"
    copy_mudata_path = tmp_path / "modified_copy.h5mu"
    original_mudata.write(str(original_mudata_path), compression="gzip")
    copied_mudata.write(str(copy_mudata_path), compression="gzip")
    return original_mudata_path, copy_mudata_path

def test_concatenate_samples_with_same_observation_ids_raises(run_component):
    """
    Test how concat handles overlapping observation IDs.
    This should raise.
    """
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
                "--input_id", "mouse;mouse2",
                "--input", input_sample1_file,
                "--input", input_sample1_file,
                "--output", "concat.h5mu",
                "--other_axis_mode", "move",
                "--output_compression", "gzip"
                ])
        assert "ValueError: Observations are not unique across samples." in \
            err.value.stdout.decode('utf-8')

@pytest.mark.parametrize("mudata_without_genome",
                          [([input_sample1_file], ["rna", "atac"])],
                          indirect=["mudata_without_genome"])
def test_concat_different_var_columns_per_sample(run_component, mudata_without_genome):
    """
    Test what happens when concatenating samples with differing auxiliary
    (like in .var) columns (present in 1 sample, absent in other).
    When concatenating the samples, all columns should be present in the
    resulting object, filling the values from samples with the missing
    column with NA.
    """
    [sample1_without_genome,] = mudata_without_genome
    run_component([
            "--input_id", "mouse;human",
            "--input", sample1_without_genome,
            "--input", input_sample2_file,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move"
            ])

    assert Path("concat.h5mu").is_file() is True
    concatenated_data = md.read("concat.h5mu")

    data_sample1 = md.read(str(sample1_without_genome))
    data_sample2 = md.read(input_sample2_file)

    assert 'genome' not in data_sample1.var_keys()
    assert concatenated_data.n_vars == data_sample1.var.index.union(data_sample2.var.index).size

    # Check if all features are present
    for mod_name in ("rna", "atac"):
        concatenated_mod = concatenated_data.mod[mod_name]
        original_var_keys = set(data_sample1.mod[mod_name].var.keys().tolist() +
                                data_sample2.mod[mod_name].var.keys().tolist())

        assert original_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in concatenated_mod.varm.keys()) | \
                            set(concatenated_mod.var.columns.tolist())
    # Value from sample1 should have NA
    assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['genome']) is True

    # Value originating from sample 2 should be the same in concatenated object.
    # Here, they are the same because we concatenate 2 samples so
    # there is only 1 unique value
    assert concatenated_data.var.loc['GL000195.1:71201-71602']['genome'] == \
            data_sample2.var.loc['GL000195.1:71201-71602']['genome']


@pytest.mark.parametrize("mudata_without_genome",
                          [([input_sample1_file, input_sample2_file], ["rna"])],
                          indirect=["mudata_without_genome"])
def test_concat_different_columns_per_modality(run_component, mudata_without_genome):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that differ between the modalities, but the difference is the same in all samples.
    """
    sample1_without_genome, sample2_without_genome = mudata_without_genome

    run_component([
            "--input_id", "mouse;human",
            "--input", sample1_without_genome,
            "--input", sample2_without_genome,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move"
            ])

    assert Path("concat.h5mu").is_file() is True
    concatenated_data = md.read("concat.h5mu")

    data_sample1 = md.read(str(sample1_without_genome))
    data_sample2 = md.read(str(sample2_without_genome))

    # Check if all features are present
    assert concatenated_data.n_vars == \
            data_sample1.var.index.union(data_sample2.var.index).size

    for mod_name in ("rna", "atac"):
        concatenated_mod = concatenated_data.mod[mod_name]
        original_var_keys = set(data_sample1.mod[mod_name].var.keys().tolist() +
                                data_sample2.mod[mod_name].var.keys().tolist())

        assert original_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in concatenated_mod.varm.keys()) | \
                            set(concatenated_mod.var.columns.tolist())

    # Check if 'interval' stays removed from modality
    assert 'genome' not in concatenated_data.mod['rna'].var.columns

    # Value from sample1 should have NA
    # ENSMUSG00000051951 is first entry from data_sample1.mod['rna'].var.head()
    assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['atac:genome']) is True
    # chr1:3094399-3095523 is the first entry from data_sample1.mod['atac'].var.head()
    assert concatenated_data.var.loc['chr1:3094399-3095523']['atac:genome'] == 'mm10'
    assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['genome'] == 'mm10'

@pytest.mark.parametrize("mudata_without_genome",
                          [([input_sample1_file], ["rna"])],
                          indirect=["mudata_without_genome"])
def test_concat_different_columns_per_modality_and_per_sample(run_component, mudata_without_genome):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that differ between the modalities and also between samples
    """

    [sample_1_without_genome, ] = mudata_without_genome
    run_component([
        "--input_id", "mouse;human",
        "--input", sample_1_without_genome,
        "--input", input_sample2_file,
        "--output", "concat.h5mu",
        "--other_axis_mode", "move"
        ])

    assert Path("concat.h5mu").is_file() == True
    concatenated_data = md.read("concat.h5mu")

    data_sample1 = md.read(str(sample_1_without_genome))
    data_sample2 = md.read(input_sample2_file)

    # Check if all features are present
    assert concatenated_data.n_vars == \
            data_sample1.var.index.union(data_sample2.var.index).size

    # Check if all features are present
    for mod_name in ("rna", "atac"):
        concatenated_mod = concatenated_data.mod[mod_name]
        original_var_keys = set(data_sample1.mod[mod_name].var.keys().tolist() +
                                data_sample2.mod[mod_name].var.keys().tolist())

        assert original_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in concatenated_mod.varm.keys()) | \
                            set(concatenated_mod.var.columns.tolist())


    # Check if 'interval' is included in RNA modality
    assert 'genome' in concatenated_data.mod['rna'].var.columns

    # Value from sample1 should have NA
    # ENSMUSG00000051951 is first entry from data_sample1.mod['rna'].var.head()
    assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['genome'])
    # chr1:3094399-3095523 is the first entry from data_sample1.mod['atac'].var.head()
    assert concatenated_data.var.loc['chr1:3094399-3095523']['genome'] == 'mm10'
    assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['genome'] == 'mm10'

@pytest.mark.parametrize("extra_column_value_sample2", [np.nan])
@pytest.mark.parametrize("extra_column_value_sample1,expected", [("bar", "bar"), (True, True), (0.1, 0.1), (np.nan, pd.NA)])
@pytest.mark.parametrize("mudata_copy_with_unique_obs",
                          [input_sample1_file],
                          indirect=["mudata_copy_with_unique_obs"])
def test_concat_remove_na(run_component, copied_mudata_with_extra_annotation_column, expected):
    """
    Test concatenation of samples where the column from one sample contains NA values
    NA values should be removed from the concatenated result
    """
    tempfile_input1, tempfile_input2 = copied_mudata_with_extra_annotation_column
    run_component([
        "--input_id", "mouse;human",
        "--input", tempfile_input1,
        "--input", tempfile_input2,
        "--output", "concat.h5mu",
        "--other_axis_mode", "move"
        ])

    assert Path("concat.h5mu").is_file() is True
    concatenated_data = md.read("concat.h5mu")
    
    if not pd.isna(expected):
        assert concatenated_data.var.loc['chr1:3094399-3095523']['test'] == expected
        assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['test'] == expected
    else:
        assert pd.isna(concatenated_data.var.loc['chr1:3094399-3095523']['test'])
        assert pd.isna(concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['test'])

    assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['test']) is True
    assert pd.isna(concatenated_data.mod['rna'].var.loc['ENSMUSG00000051951']['test']) is True

@pytest.mark.parametrize("extra_column_annotation_matrix", ["obs"])
@pytest.mark.parametrize("extra_column_value_sample1,extra_column_value_sample2,expected", [(1, "1", pd.CategoricalDtype(categories=['1.0', '1']))])
@pytest.mark.parametrize("mudata_copy_with_unique_obs",
                          [input_sample1_file],
                          indirect=["mudata_copy_with_unique_obs"])
def test_concat_dtypes(run_component, copied_mudata_with_extra_annotation_column, expected):
    """
    Test joining column with different dtypes to make sure that they are writable.
    The default path is to convert all non-na values to strings and wrap the column into a categorical dtype.
    """
    tempfile_input1, tempfile_input2 = copied_mudata_with_extra_annotation_column
    run_component([
        "--input_id", "mouse;human",
        "--input", tempfile_input1,
        "--input", tempfile_input2,
        "--output", "concat.h5mu",
        "--other_axis_mode", "move"
        ])
    concatenated_data = md.read("concat.h5mu")
    concatenated_data.mod['atac'].obs['test'].dtype == expected

@pytest.mark.parametrize("extra_column_annotation_matrix", ["var"])
@pytest.mark.parametrize("extra_column_value_sample1,extra_column_value_sample2", [("2", "1")])
@pytest.mark.parametrize("mudata_copy_with_unique_obs",
                          [input_sample1_file],
                          indirect=["mudata_copy_with_unique_obs"])
def test_resolve_annotation_conflict_missing_column(run_component, copied_mudata_with_extra_annotation_column, make_obs_names_unique, tmp_path):
    """
    Test using mode 'move' and resolving a conflict in metadata between the samples,
    but the metadata column is missing in one of the samples.
    """
    tempfile_input1, tempfile_input2 = copied_mudata_with_extra_annotation_column
    original_data = md.read_h5mu(input_sample1_file)
    make_obs_names_unique(original_data)
    original_data_path = tmp_path / f"{uuid.uuid4().hex}.h5mu"
    original_data.write_h5mu(original_data_path)
    run_component([
        "--input_id", "mouse;human;sample_without_column",
        "--input", tempfile_input1,
        "--input", tempfile_input2,
        "--input", original_data_path,
        "--output", "concat.h5mu",
        "--other_axis_mode", "move"
        ])
    concatenated_data = md.read("concat.h5mu")
    assert 'test' not in concatenated_data.mod['atac'].var.columns
    assert 'test' not in concatenated_data.mod['atac'].obs.columns
    assert 'conflict_test' in concatenated_data.mod['atac'].varm

def test_mode_move(run_component, tmp_path):
    tempfile_input1 = tmp_path / "input1.h5mu"
    tempfile_input2 = tmp_path / "input2.h5mu"
    input1 = md.read(input_sample1_file)
    input2 = md.read(input_sample2_file)
    # Create conflict
    var_index_values = input1.var.index.to_numpy()
    var_index_values[0] = input2.var.index[0]
    input1.var.index = var_index_values
    var_index_values = input1.mod['rna'].var.index.to_numpy()
    var_index_values[0] = input2.mod['rna'].var.index.to_numpy()[0]
    input1.mod['rna'].var.index = var_index_values
    input1.write(tempfile_input1.name)
    input2.write(tempfile_input2.name)
    run_component([
        "--input_id", "mouse;human",
        "--input", tempfile_input1.name,
        "--input", tempfile_input2.name,
        "--output", "concat.h5mu",
        "--other_axis_mode", "move"
        ])
    assert Path("concat.h5mu").is_file() is True
    concatenated_data = md.read("concat.h5mu")

    # Check if observations from all of the samples are present
    assert (concatenated_data.n_obs ==  input1.n_obs + input2.n_obs)

    # Check if all modalities are present
    sample1_mods, sample2_mods = set(input1.mod.keys()), set(input2.mod.keys())
    concatentated_mods = set(concatenated_data.mod.keys())
    assert (sample1_mods | sample2_mods) == concatentated_mods

    varm_check = {
        "rna": ({"conflict_gene_symbol": ("mouse", "human"),
                 "conflict_genome": ("mouse", "human")}),
        "atac": {}
    }

    # Check if all features are present
    for mod_name in ("rna", "atac"):
        concatenated_mod = concatenated_data.mod[mod_name]
        original_var_keys = set(input1.mod[mod_name].var.keys().tolist() +
                                    input2.mod[mod_name].var.keys().tolist())

        assert original_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in concatenated_mod.varm.keys()) | \
                            set(concatenated_mod.var.columns.tolist())

        varm_expected = varm_check[mod_name]
        assert list(concatenated_mod.varm.keys()) == list(varm_expected.keys())
        for varm_key, expected_columns in varm_expected.items():
            assert tuple(concatenated_mod.varm[varm_key].columns) == expected_columns
        if not varm_expected:
            assert concatenated_mod.varm == {}
        assert concatenated_mod.obsm == {}

def test_concat_invalid_h5_error_includes_path(run_component, tmp_path):
    empty_file = tmp_path / "empty.h5mu"
    empty_file.touch()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
                "--input_id", "mouse;empty",
                "--input", input_sample1_file,
                "--input", empty_file,
                "--output", "concat.h5mu",
                "--other_axis_mode", "move"
                ])
        assert re.search(rf"OSError: Failed to load .*{str(empty_file)}\. Is it a valid h5 file?",
            err.value.stdout.decode('utf-8'))
        

@pytest.mark.parametrize("mudata_without_genome",
                          [([input_sample1_file], ["rna", "atac"])],
                          indirect=["mudata_without_genome"])
def test_concat_var_obs_names_order(run_component, mudata_without_genome, 
                                    anndata_to_sparse_dataframe):
    """
    Test that the var_names and obs_names are still linked to the correct count data.
    """
    [sample1_without_genome,] = mudata_without_genome
    run_component([
            "--input_id", "mouse;human",
            "--input", sample1_without_genome,
            "--input", input_sample2_file,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move"
            ])
    assert Path("concat.h5mu").is_file() is True
    for sample_name, sample_path in {"mouse": sample1_without_genome, 
                                     "human": input_sample2_file}.items():
        for mod_name in ["rna", "atac"]:
            data_sample = md.read_h5ad(sample_path, mod=mod_name)
            processed_data = md.read_h5ad("concat.h5mu", mod=mod_name)
            muon.pp.filter_obs(processed_data, 'sample_id', lambda x: x == sample_name)
            muon.pp.filter_var(processed_data, data_sample.var_names)
            data_sample_to_test = anndata_to_sparse_dataframe(data_sample)
            processed_data_to_test = anndata_to_sparse_dataframe(processed_data)
            data_sample_to_test = data_sample_to_test.reindex_like(processed_data_to_test)
            pd.testing.assert_index_equal(data_sample_to_test.columns, processed_data_to_test.columns)
            pd.testing.assert_index_equal(data_sample_to_test.index, processed_data_to_test.index)
            for (_, col1), (_, col2) in zip(data_sample_to_test.items(), processed_data_to_test.items()):
                pd._testing.assert_sp_array_equal(col1.array, col2.array)
            

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "-v"]))
