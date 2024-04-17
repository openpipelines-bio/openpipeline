# import os
from mudata import read_h5mu
# import pytest
# import sys
from openpipelinetestutils.asserters import assert_annotation_objects_equal
from openpipelinetestutils.utils import remove_annotation_column


##VIASH START
par = {
  "input": "/Users/jakubmajercik/Data_Intuitive/openpipelines-bio/openpipeline/work/8b/998dbb9cc41bb91c78eebf17596716/mouse_test.publish.output.h5mu",
  "og_input": "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"
}

meta = {
    "resources_dir": "resources_test/concat_test_data",
}
##VIASH END
    
    
# def test_fields_in_output():
input_mudata = read_h5mu(par["og_input"])
output_mudata = read_h5mu(par["input"])

assert input_mudata.n_mod == output_mudata.n_mod, "Number of modalities should be the same"
assert input_mudata.mod.keys() == output_mudata.mod.keys(), "Modalities should be the same"
assert list(output_mudata.mod.keys()) == ["rna", "atac"], "Modalities should be rna and atac"

obs_cols_to_remove = []
for top_n_vars in ("50", "100", "200", "500"):
    obs_cols_to_remove.append(f"pct_of_counts_in_top_{top_n_vars}_vars")

obs_cols_to_remove.extend(['total_counts', 'num_nonzero_vars'])
var_cols_to_remove = ['obs_mean', 'total_counts', 'num_nonzero_obs', 'pct_dropout']

assert set(obs_cols_to_remove).issubset(set(output_mudata.mod["rna"].obs.columns.to_list()))
assert set(var_cols_to_remove).issubset(set(output_mudata.mod["rna"].var.columns.to_list()))

initial_mudata = remove_annotation_column(output_mudata, obs_cols_to_remove, axis="obs", modality_name="rna")
initial_mudata = remove_annotation_column(initial_mudata, var_cols_to_remove, axis="var", modality_name="rna")

assert_annotation_objects_equal(input_mudata, initial_mudata)

    
# if __name__ == "__main__":
#     sys.exit(pytest.main([__file__]))