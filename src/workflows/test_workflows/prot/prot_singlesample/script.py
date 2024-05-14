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

assert_annotation_objects_equal(input_mudata.mod["rna"], output_mudata.mod["rna"])

assert 'filter_with_counts' in output_mudata.mod["prot"].var.keys()