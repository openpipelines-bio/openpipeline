import mudata as mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
  "input": "output_test/split_modalities/foo_types.csv",
  "orig_input": "output_test/split_modalities/h5mu"
}

meta = {
    "resources_dir": "resources_test/pbmc_1k_protein_v3"
}

##VIASH END

print ("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"
assert_annotation_objects_equal(input.mod["atac"], output.mod["atac"])




print("Test successful!", flush=True)