from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
    "input": "input.h5mu",
    "input": "input_og.h5mu",
    "is_corrected": True
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input_og"])
output_mudata = read_h5mu(par["input"])

assert input_mudata.mod.keys() == output_mudata.mod.keys(), "Input and output should have the same modalities."

for modality,input_adata,output_adata in zip(input_mudata.mod.keys(),
                                             input_mudata.mod.values(),
                                             output_mudata.mod.values()):
    assert input_adata.n_obs >= output_adata.n_obs, "Output should have less or equal number of observations than input."
    assert input_adata.n_vars == output_adata.n_vars, "Output should have the same number of variables as input."
    if modality != "rna":
        assert_annotation_objects_equal(input_adata,
                                        output_adata)

if par["is_corrected"]:
    assert "cellbender_corrected" in output_mudata.mod["rna"].layers