from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
    "input": "input.h5mu",
    "input": "input_og.h5mu",
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input_og"])
output_mudata = read_h5mu(par["input"])

assert input_mudata.mod.keys() == output_mudata.mod.keys(), "Input and output should have the same modalities."

assert_annotation_objects_equal(input_mudata.mod["atac"], output_mudata.mod["atac"])

assert set(['num_nonzero_vars',
            'total_counts',
            'pct_of_counts_in_top_50_vars',
            'pct_of_counts_in_top_100_vars',
            'pct_of_counts_in_top_200_vars',
            'pct_of_counts_in_top_500_vars']).issubset(output_mudata.mod["rna"].obs.keys())
assert set(['filter_with_hvg',
            'obs_mean',
            'total_counts',
            'num_nonzero_obs',
            'pct_dropout']).issubset(output_mudata.mod["rna"].var.keys())
assert 'log_normalized' in output_mudata.mod["rna"].layers.keys()
assert 'log1p' in output_mudata.mod["rna"].uns.keys()