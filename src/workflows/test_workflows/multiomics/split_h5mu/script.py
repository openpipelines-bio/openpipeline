import csv
import os
import mudata as mu
# from openpipeline_testutils.asserters import assert_annotation_objects_equal


##VIASH START
par = {
    "input": "/home/di/code/openpipelines-multisample/work/82/485f6cde9bab8a53f15726b59de4c1/_viash_par/input_1/samples.csv",
    "samples_dir": "/home/di/code/openpipelines-multisample/work/82/485f6cde9bab8a53f15726b59de4c1/_viash_par/samples_dir_1/samples",
    "orig_input": "/home/di/code/openpipelines-multisample/work/82/485f6cde9bab8a53f15726b59de4c1/_viash_par/orig_input_1/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END


print("Loading data", flush=True)
with open(par["input"], "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    data = list(reader)

input_mu = mu.read_h5mu(par["orig_input"])

input_samples = input_mu.mod["rna"].obs["sample_id"].unique()
num_samples = input_samples.shape[0]
num_files = len(os.listdir(par["samples_dir"]))

# Check if the number of files is equal to the number of lines in the csv
assert num_samples == num_files, f"Expected {num_samples} files, but found {num_files}."
assert (
    input_mu.n_mod == num_samples
), f"Expected {num_samples} samples in {par['orig_input']} got {num_files} samples."


# Check if the files exist and if the sample name is in the file name
for i, row in enumerate(data):
    if i == 0:
        continue
    # Check if the files exist and if the sample name is in the file name
    assert row[0] in row[1], f"Expected {row[0]} to be in {row[1]}."
    samples_fp = os.path.join(par["samples_dir"], row[1])
    assert os.path.exists(samples_fp), f"Expected {row[1]} to exist."
    # Check sample is correct in the h5mu file
    mod_mu = mu.read_h5mu(samples_fp)
    found_samples = mod_mu["rna"].obs["sample_id"].unique()
    assert found_samples.shape[0] == 1, f"Expected 1 sample ID in {row[1]}."
    assert (
        found_samples[0] == row[0]
    ), f"Expected {row[0]} to be the sample ID for {row[1]}"
    assert (
        row[0] in input_samples
    ), f"Expected {row[0]} to be a mod in {par['orig_input']}."

print("Test successful!", flush=True)
