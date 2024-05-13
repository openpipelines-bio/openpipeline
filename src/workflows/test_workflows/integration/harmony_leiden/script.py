from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
    "input": "input.h5mu",
    "input": "input_og.h5mu",
    "leiden_resolution": [1]
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input_og"])
output_mudata = read_h5mu(par["input"])

assert input_mudata.mod.keys() == output_mudata.mod.keys(), "Input and output should have the same modalities."

assert "bbknn_integration_neighbors" in output_mudata.mod["rna"].uns, "Output should contain key 'harmonypy_integration_neighbors'."
assert "X_leiden_bbknn_umap" in output_mudata.mod["rna"].obsm, "Output should contain key 'X_leiden_harmony_umap'."
assert set(['harmonypy_integration_connectivities', 'harmonypy_integration_distances']).issubset(output_mudata.mod["rna"].obsp.keys()), "Output should contain keys 'bbknn_integration_connectivities' and 'bbknn_integration_distances'."

# if not par["leiden_resolution"]:
#     assert "bbknn_integration_leiden_1.0" in output_mudata.mod["rna"].obs, "Output should contain key 'bbknn_integration_leiden_1.0'."

assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])