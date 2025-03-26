import subprocess
from os import path
import mudata as mu
import numpy as np

## VIASH START
meta = {
    "executable": "target/docker/convert/from_bdrhap_to_h5mu/from_bdrhap_to_h5mu",
    "resources_dir": "resources_test/bdrhap_5kjrt/processed/output_raw/",
}
## VIASH END

input = meta["resources_dir"] + "/sample.h5mu"
output = "output1.h5mu"

cmd_pars = [
    meta["executable"],
    "--input",
    input,
    "--output",
    output,
    "--id",
    "foo",
    "--output_compression",
    "gzip",
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

print(">> Check if output exists", flush=True)
assert path.exists(output), "No output was created."


print(">> Check contents of output", flush=True)
data = mu.read_h5mu(output)
rna_adata = data.mod["rna"]
prot_adata = data.mod["prot"]

# check whether correct feature types are detected
assert np.array_equal(
    rna_adata.var["feature_type"].unique(), ["Gene Expression"]
), "RNA expression should only contain Gene Expression vars."
assert np.array_equal(
    rna_adata.var["reference_file"].unique(), ["reference_bd_rhapsody.tar.gz"]
), "Wrong reference file detected for Gene Expression vars."
assert "ADAMTSL4" in rna_adata.var_names, 'RNA modality should contain gene "ADAMTS4".'
assert np.array_equal(
    rna_adata.obs["library_id"].unique(), ["12ABC & 12SMK & 12WTA"]
), "Gene Expression .obs library_id should equal '12ABC & 12WTA."
assert (
    "sample_tag" in rna_adata.obs.keys()
), "RNA modality should contain column 'sample_id'."
assert (
    "sample_id" in rna_adata.obs.keys()
), "RNA modality should contain column 'sample_name'."

assert np.array_equal(
    prot_adata.var["feature_type"].unique(), ["Antibody Capture"]
), "RNA expression should only contain Antibody Capture vars."
assert np.array_equal(
    prot_adata.var["reference_file"].unique(), ["BDAbSeq_ImmuneDiscoveryPanel.fasta"]
), "Wrong reference file detected for Antibody Capture vars."
assert (
    "CD279:EH12-1|PDCD1|AHS0014|pAbO" in prot_adata.var_names
), 'Protein modality should contain protein "CD279:EH12-1|PDCD1|AHS0014|pAbO".'
assert np.array_equal(
    prot_adata.obs["library_id"].unique(), ["12ABC & 12SMK & 12WTA"]
), "Antibody Capture .obs library_id should equal '12ABC & 12WTA."
assert (
    "sample_tag" in prot_adata.obs.keys()
), "Protein modality should contain column 'sample_id'."
assert (
    "sample_id" in prot_adata.obs.keys()
), "Protein modality should contain column 'sample_name'."

# check whether gene was found
assert "PDE4DIP" in data.var_names, 'Output should contain gex column "PDE4DIP".'
assert (
    "CD279:EH12-1|PDCD1|AHS0014|pAbO" in data.var_names
), 'Output should contain abc column "CD279:EH12-1|PDCD1|AHS0014|pAbO".'

print("> Test successful", flush=True)
