import csv
import os
import mudata as mu
from openpipeline_testutils.asserters import assert_annotation_objects_equal


##VIASH START
par = {
    "input": "output_test/split_modalities/foo_types.csv",
    "mod_dir": "output_test/split_modalities/h5mu",
    "orig_input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3.h5mu",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END


print("Loading data", flush=True)
with open(par["input"], "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    data = list(reader)

input_mu = mu.read_h5mu(par["orig_input"])

num_mod = len(data)
num_files = len(os.listdir(par["mod_dir"]))

# Check if the number of files is equal to the number of lines in the csv
assert num_mod == num_files, f"Expected {num_mod} files, but found {num_files}."
assert input_mu.n_mod == num_mod, (
    f"Expected {num_mod} modalities in {par['orig_input']} got {input_mu.n_mod} modalities."
)

output_mods = {
    csv_entry["name"]: mu.read_h5mu(os.path.join(par["mod_dir"], csv_entry["filename"]))
    for csv_entry in data
}

# Check if the files exist and if the modality name is in the file name
for csv_item in data:
    mod_name, mudata_file_name = csv_item["name"], csv_item["filename"]
    # Check if the files exist and if the modality name is in the file name
    assert mod_name in mudata_file_name, (
        f"Expected {mod_name} to be in {mudata_file_name}."
    )
    mudata_path = os.path.join(par["mod_dir"], mudata_file_name)
    assert os.path.exists(mudata_path), f"Expected {mudata_file_name} to exist."
    # Check modality is correct in the h5mu file
    mudata_object = mu.read_h5mu(mudata_path)
    assert mudata_object.n_mod == 1, f"Expected 1 modality in {mudata_file_name}."
    assert mod_name in mudata_object.mod.keys(), (
        f"Expected {mod_name} to be the mod in {mudata_file_name}."
    )
    assert mod_name in input_mu.mod.keys(), (
        f"Expected {mod_name} to be a mod in {par['orig_input']}."
    )
    # Check if extracted modalities are equal to the original modalities
    assert_annotation_objects_equal(mudata_object[mod_name], input_mu.mod[mod_name])

print("Test successful!", flush=True)
