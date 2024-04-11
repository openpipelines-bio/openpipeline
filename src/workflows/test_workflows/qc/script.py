import os
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal


##VIASH START
par = {
  "input": "/Users/jakubmajercik/Data_Intuitive/openpipelines-bio/openpipeline/work/8b/998dbb9cc41bb91c78eebf17596716/mouse_test.publish.output.h5mu",
  "og_input": "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"
}

meta = {
    "resources_dir": "resources_test/concat_test_data",
}
##VIASH END

print("Here")
assert True
# input_mudata = read_h5mu(par["og_input"])
# output_mudata = read_h5mu(par["input"])

# for top_n_vars in ("50", "100", "200", "500"):
#     assert f"pct_of_counts_in_top_{top_n_vars}_vars" in output_mudata.mod['rna'].obs