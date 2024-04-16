import mudata as mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
  "input": "resources_test/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu",
  "orig_input": "test.h5mu"
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

# Check vdj_t modality
assert_annotation_objects_equal(input.mod["vdj_t"], output.mod["vdj_t"])

# Check prot modality
assert_annotation_objects_equal(input.mod["prot"], output.mod["prot"])

# Check rna modality
assert_annotation_objects_equal(input.mod["rna"], output.mod["rna"])

print("Test successful!", flush=True)