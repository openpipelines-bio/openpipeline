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

assert "totalvi_integration_leiden_1.0" in output_mudata.mod["rna"].obs, "Output should contain key 'totalvi_integration_leiden_1.0'."
assert "totalvi_integration_neighbors" in output_mudata.mod["rna"].uns, "Output should contain key 'totalvi_integration_neighbors'."
assert set(['X_totalvi_umap', 'X_integrated_totalvi', 'X_totalvi']).issubset(output_mudata.mod["rna"].obsm.keys())
assert set(['totalvi_integration_connectivities', 'totalvi_integration_distances']).issubset(output_mudata.mod["rna"].obsp.keys()), "Output should contain keys 'scvi_integration_connectivities' and 'scvi_integration_distances'."

assert "X_totalvi" in output_mudata.mod["prot"].obsm.keys()