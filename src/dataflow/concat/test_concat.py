import mudata as md
import subprocess
from pathlib import Path
import pandas as pd
from tempfile import NamedTemporaryFile
import numpy as np
import pytest

## VIASH START
meta = {
    'executable': './target/docker/dataflow/concat/concat',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2
}
## VIASH END

meta['cpus'] = 1 if not meta['cpus'] else meta['cpus']

resources_dir = meta["resources_dir"]
# Note: the .var for these samples have no overlap, so there are no conflicting annotations
# for the features that need to be handled by the concat component.
# The tests below that specifically test the concatenation of conflicting data need to introduce
# the conflict.
input_sample1_file = f"{resources_dir}/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"
input_sample2_file = f"{resources_dir}/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"

def test_concatenate_samples_with_same_observation_ids_raises(run_component):
    """
    Test how concat handles overlapping observation IDs.
    This should raise.
    """
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
                "--input_id", "mouse,mouse2",
                "--input", input_sample1_file,
                "--input", input_sample1_file,
                "--output", "concat.h5mu",
                "--other_axis_mode", "move",
                "---cpus", str(meta["cpus"])
                ])
    err.match(r"ValueError: Observations are not unique across samples\.")

def test_concat_different_var_columns_per_sample(run_component):
    """
    Test what happens when concatenating samples with differing auxiliary
    (like in .var) columns (present in 1 sample, absent in other).
    When concatenating the samples, all columns should be present in the
    resulting object, filling the values from samples with the missing
    column with NA.
    """

    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile:
        data_sample1 = md.read(input_sample1_file)
        data_sample1.mod['rna'].var.drop('genome', axis=1, inplace=True)
        data_sample1.mod['atac'].var.drop('genome', axis=1, inplace=True)
        data_sample1.var.drop('genome', axis=1, inplace=True)
        data_sample1.write(tempfile.name, compression="gzip")
        tempfile.flush()
        run_component([
                "--input_id", "mouse,human",
                "--input", tempfile.name,
                "--input", input_sample2_file,
                "--output", "concat.h5mu",
                "--other_axis_mode", "move",
                "---cpus", str(meta["cpus"])
                ])

        assert Path("concat.h5mu").is_file() is True
        concatenated_data = md.read("concat.h5mu")

        data_sample1 = md.read(tempfile.name)
        data_sample2 = md.read(input_sample2_file)

        assert 'genome' not in data_sample1.var_keys()
        assert concatenated_data.n_vars == data_sample1.var.index.union(data_sample2.var.index).size
        # Check if all features are present
        rna = concatenated_data.mod['rna']
        atac = concatenated_data.mod['atac']
        original_rna_var_keys = set(data_sample1.mod['rna'].var.keys().tolist() +
                                    data_sample2.mod['rna'].var.keys().tolist())
        original_atac_var_keys = set(data_sample1.mod['atac'].var.keys().tolist() + 
                                    data_sample2.mod['atac'].var.keys().tolist())
        assert original_rna_var_keys == \
                            set(column_name.removeprefix('conflict_') 
                                for column_name in rna.varm.keys()) | \
                            set(rna.var.columns.tolist())
        assert original_atac_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in atac.varm.keys()) | \
                            set(atac.var.columns.tolist())
        # Value from sample1 should have NA
        assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['genome']) is True

        # Value originating from sample 2 should be the same in concatenated object.
        # Here, they are the same because we concatenate 2 samples so
        # there is only 1 unique value
        assert concatenated_data.var.loc['GL000195.1:71201-71602']['genome'] == \
               data_sample2.var.loc['GL000195.1:71201-71602']['genome']

def test_concat_different_columns_per_modality(run_component):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that differ between the modalities, but the difference is the same in all samples.
    """
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_sample1,\
            NamedTemporaryFile('w', suffix=".h5mu") as tempfile_sample2:
        data_sample1 = md.read(input_sample1_file)
        data_sample1.mod['rna'].var.drop('genome', axis=1, inplace=True)
        data_sample2 = md.read(input_sample2_file)
        data_sample2.mod['rna'].var.drop('genome', axis=1, inplace=True)
        # Note: when a column is present in the feature annotation for one
        # modality (MuData.mod[...].var) but not for another other, then the global
        # var (MuData.var) gets a new column 'mod_name:column_name'
        # (here atac:genome) next to the old 'column_name' (here just 'genome')
        data_sample1.update_var()
        data_sample2.update_var()
        data_sample1.var.drop('genome', axis=1, inplace=True)
        data_sample2.var.drop('genome', axis=1, inplace=True)
        data_sample1.write(tempfile_sample1.name, compression="gzip")
        data_sample2.write(tempfile_sample2.name, compression="gzip")

        run_component([
                "--input_id", "mouse,human",
                "--input", tempfile_sample1.name,
                "--input", tempfile_sample2.name,
                "--output", "concat.h5mu",
                "--other_axis_mode", "move",
                "---cpus", str(meta["cpus"])
                ])

        assert Path("concat.h5mu").is_file() is True
        concatenated_data = md.read("concat.h5mu")

        data_sample1 = md.read(tempfile_sample1.name)
        data_sample2 = md.read(tempfile_sample2.name)

        # Check if all features are present
        assert concatenated_data.n_vars == \
               data_sample1.var.index.union(data_sample2.var.index).size

        # Check if all features are present
        rna = concatenated_data.mod['rna']
        atac = concatenated_data.mod['atac']
        original_rna_var_keys = set(data_sample1.mod['rna'].var.keys().tolist() +
                                    data_sample2.mod['rna'].var.keys().tolist())
        original_atac_var_keys = set(data_sample1.mod['atac'].var.keys().tolist() + 
                                    data_sample2.mod['atac'].var.keys().tolist())
        assert original_rna_var_keys == \
                            set(column_name.removeprefix('conflict_') 
                                for column_name in rna.varm.keys()) | \
                            set(rna.var.columns.tolist())
        assert original_atac_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in atac.varm.keys()) | \
                            set(atac.var.columns.tolist())

        # Check if 'interval' stays removed from modality
        assert 'genome' not in concatenated_data.mod['rna'].var.columns

        # Value from sample1 should have NA
        # ENSMUSG00000051951 is first entry from data_sample1.mod['rna'].var.head()
        assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['atac:genome']) is True
        # chr1:3094399-3095523 is the first entry from data_sample1.mod['atac'].var.head()
        assert concatenated_data.var.loc['chr1:3094399-3095523']['atac:genome'] == 'mm10'
        assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['genome'] == 'mm10'

def test_concat_different_columns_per_modality_and_per_sample(run_component):
    """
    Test what happens when concatenating samples that have auxiliary columns
    that differ between the modalities and also between samples
    """
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_sample1:
        data_sample1 = md.read(input_sample1_file)
        data_sample1.mod['rna'].var.drop('genome', axis=1, inplace=True)
        # Note: when a column is present in the feature annotation for one
        # modality (MuData.mod[...].var) but not for another other, then the global
        # var (MuData.var) gets a new column 'mod_name:column_name'
        # (here atac:genome) next to the old 'column_name' (here just 'genome')
        data_sample1.update_var()
        data_sample1.var.drop('genome', axis=1, inplace=True)
        data_sample1.write(tempfile_sample1.name, compression="gzip")

        run_component([
            "--input_id", "mouse,human",
            "--input", tempfile_sample1.name,
            "--input", input_sample2_file,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move",
            "---cpus", str(meta["cpus"])
            ])

        assert Path("concat.h5mu").is_file() == True
        concatenated_data = md.read("concat.h5mu")

        data_sample1 = md.read(tempfile_sample1.name)
        data_sample2 = md.read(input_sample2_file)

        # Check if all features are present
        assert concatenated_data.n_vars == \
               data_sample1.var.index.union(data_sample2.var.index).size

        # Check if all features are present
        rna = concatenated_data.mod['rna']
        atac = concatenated_data.mod['atac']
        original_rna_var_keys = set(data_sample1.mod['rna'].var.keys().tolist() +
                                    data_sample2.mod['rna'].var.keys().tolist())
        original_atac_var_keys = set(data_sample1.mod['atac'].var.keys().tolist() + 
                                    data_sample2.mod['atac'].var.keys().tolist())
        assert original_rna_var_keys == \
                            set(column_name.removeprefix('conflict_') 
                                for column_name in rna.varm.keys()) | \
                            set(rna.var.columns.tolist())
        assert original_atac_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in atac.varm.keys()) | \
                            set(atac.var.columns.tolist())

        # Check if 'interval' is included in RNA modality
        assert 'genome' in concatenated_data.mod['rna'].var.columns

        # Value from sample1 should have NA
        # ENSMUSG00000051951 is first entry from data_sample1.mod['rna'].var.head()
        assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['genome'])
        # chr1:3094399-3095523 is the first entry from data_sample1.mod['atac'].var.head()
        assert concatenated_data.var.loc['chr1:3094399-3095523']['genome'] == 'mm10'
        assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['genome'] == 'mm10'


def test_concat_remove_na(run_component):
    """
    Test concatenation of samples where the column from one sample contains NA values
    NA values should be removed from the concatenated result
    """
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input1,\
            NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input2:
        input1 = md.read(input_sample1_file)
        input2 = input1.copy()
        input2.obs.index = "make_unique_" + input2.obs.index
        for mod_data in input2.mod.values():
            mod_data.obs.index = "make_unique_" + mod_data.obs.index            
        input1.mod['rna'].var['test'] = np.nan
        input1.mod['atac'].var['test'] = "bar"
        input1.update_var()
        input2.mod['rna'].var['test'] = np.nan
        input2.mod['atac'].var['test'] = np.nan
        input2.update_var()
        input2.update_obs()
        input1.write(tempfile_input1.name, compression="gzip")
        input2.write(tempfile_input2.name, compression="gzip")

        run_component([
            "--input_id", "mouse,human",
            "--input", tempfile_input1.name,
            "--input", tempfile_input2.name,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move",
            "---cpus", str(meta["cpus"])
            ])

        assert Path("concat.h5mu").is_file() is True
        concatenated_data = md.read("concat.h5mu")

        assert concatenated_data.var.loc['chr1:3094399-3095523']['test'] == 'bar'
        assert  concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['test'] == 'bar'

        assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['test']) is True
        assert pd.isna(concatenated_data.mod['rna'].var.loc['ENSMUSG00000051951']['test']) is True


def test_concat_only_bool(run_component):
    """
    Test concatenation of samples where the column from one sample contains NA values
    NA values should be removed from the concatenated result
    """
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input1,\
            NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input2:
        input1 = md.read(input_sample1_file)
        input2 = input1.copy()
        input2.obs.index = "make_unique_" + input2.obs.index
        for mod_data in input2.mod.values():
            mod_data.obs.index = "make_unique_" + mod_data.obs.index
        input1.mod['rna'].var['test'] = np.nan
        input1.mod['atac'].var['test'] = True
        input1.update()
        input2.mod['rna'].var['test'] = np.nan
        input2.mod['atac'].var['test'] = np.nan
        input2.update()
        input1.write(tempfile_input1.name, compression="gzip")
        input2.write(tempfile_input2.name, compression="gzip")

        run_component([
            "--input_id", "mouse,human",
            "--input", tempfile_input1.name,
            "--input", tempfile_input2.name,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move",
            "---cpus", str(meta["cpus"])
            ])

        assert Path("concat.h5mu").is_file() is True
        concatenated_data = md.read("concat.h5mu")

        assert concatenated_data.var.loc['chr1:3094399-3095523']['test'] == 'True'
        assert concatenated_data.mod['atac'].var.loc['chr1:3094399-3095523']['test'] == 'True'

        assert pd.isna(concatenated_data.var.loc['ENSMUSG00000051951']['test']) is True
        assert pd.isna(concatenated_data.mod['rna'].var.loc['ENSMUSG00000051951']['test']) is True

def test_mode_move(run_component):
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input1,\
            NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input2:
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
            "--input_id", "mouse,human",
            "--input", tempfile_input1.name,
            "--input", tempfile_input2.name,
            "--output", "concat.h5mu",
            "--other_axis_mode", "move",
            "---cpus", str(meta["cpus"])]
            )
        assert Path("concat.h5mu").is_file() is True
        concatenated_data = md.read("concat.h5mu")

        data_sample1 = md.read(tempfile_input1.name)
        data_sample2 = md.read(tempfile_input2.name)

        # Check if observations from all of the samples are present
        assert (concatenated_data.n_obs ==  data_sample1.n_obs + data_sample2.n_obs)

        # Check if all features are present
        rna = concatenated_data.mod['rna']
        atac = concatenated_data.mod['atac']
        original_rna_var_keys = set(data_sample1.mod['rna'].var.keys().tolist() +
                                    data_sample2.mod['rna'].var.keys().tolist())
        original_atac_var_keys = set(data_sample1.mod['atac'].var.keys().tolist() + 
                                    data_sample2.mod['atac'].var.keys().tolist())
        assert original_rna_var_keys == \
                            set(column_name.removeprefix('conflict_') 
                                for column_name in rna.varm.keys()) | \
                            set(rna.var.columns.tolist())
        assert original_atac_var_keys == \
                            set(column_name.removeprefix('conflict_')
                                for column_name in atac.varm.keys()) | \
                            set(atac.var.columns.tolist())

        # Check if all modalities are present
        sample1_mods, sample2_mods = set(data_sample1.mod.keys()), set(data_sample2.mod.keys())
        concatentated_mods = set(concatenated_data.mod.keys())
        assert (sample1_mods | sample2_mods) == concatentated_mods

        # Check if conflicting columns in .varm
        assert  list(rna.varm.keys()) == ["conflict_gene_symbol", "conflict_genome"]
        assert  rna.varm['conflict_gene_symbol'].columns.tolist() == ["mouse", "human"]
        assert rna.varm['conflict_genome'].columns.tolist() == ["mouse", "human"]
        assert dict(atac.varm) == {}
        assert dict(rna.obsm) == {}
        assert dict(atac.obsm) == {}


if __name__ == '__main__':
    pytest.main([__file__], plugins=["viashpy"])
