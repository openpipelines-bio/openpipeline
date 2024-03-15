import csv
import os


##VIASH START
par = {
  "input": "output_test/split_modalities/foo_types.csv",
  "mod_dir": "output_test/split_modalities/h5mu"
}

meta = {
    "resources_dir": "resources_test/pbmc_1k_protein_v3"
}

##VIASH END


print("Loading data", flush=True)
with open(par["input"], "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    data = list(reader)

csv_lines = len(data) -1
num_files = len(os.listdir(par["mod_dir"]))

# Check if the number of files is equal to the number of lines in the csv
assert csv_lines == num_files, f"Expected {csv_lines} files, but found {num_files}."

# Check if the files exist and if the modality name is in the file name
for i, row in enumerate(data):
    if i == 0:
        continue
    assert row[0] in row[1], f"Expected {row[0]} to be in {row[1]}."
    mod_fp = os.path.join(par["mod_dir"], row[1])
    assert os.path.exists(mod_fp), f"Expected {row[1]} to exist."
    # Check modality is correct in the h5mu file
    

print("Test successful!", flush=True)
