import subprocess
from os import path
import mudata as mu
import numpy as np
import awkward as ak

## VIASH START
meta = {
    "executable": "target/docker/convert/from_bdrhap_to_h5mu/from_bdrhap_to_h5mu",
    "resources_dir": "resources_test/bdrhap_5kjrt/processed/output_raw/",
}
## VIASH END

input = meta["resources_dir"] + "/sample.cellismo"
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
assert np.array_equal(rna_adata.var["feature_type"].unique(), ["Gene Expression"]), (
    "RNA expression should only contain Gene Expression vars."
)
assert np.array_equal(
    rna_adata.var["reference_file"].unique(), ["reference_bd_rhapsody.tar.gz"]
), "Wrong reference file detected for Gene Expression vars."
assert "ADAMTSL4" in rna_adata.var_names, 'RNA modality should contain gene "ADAMTS4".'
assert np.array_equal(
    rna_adata.obs["library_id"].unique(), ["12ABC & 12SMK & 12WTA"]
), "Gene Expression .obs library_id should equal '12ABC & 12WTA."
assert "sample_tag" in rna_adata.obs.keys(), (
    "RNA modality should contain column 'sample_id'."
)
assert "sample_id" in rna_adata.obs.keys(), (
    "RNA modality should contain column 'sample_name'."
)

assert np.array_equal(prot_adata.var["feature_type"].unique(), ["Antibody Capture"]), (
    "RNA expression should only contain Antibody Capture vars."
)
assert np.array_equal(
    prot_adata.var["reference_file"].unique(), ["BDAbSeq_ImmuneDiscoveryPanel.fasta"]
), "Wrong reference file detected for Antibody Capture vars."
assert "CD279:EH12-1|PDCD1|AHS0014|pAbO" in prot_adata.var_names, (
    'Protein modality should contain protein "CD279:EH12-1|PDCD1|AHS0014|pAbO".'
)
assert np.array_equal(
    prot_adata.obs["library_id"].unique(), ["12ABC & 12SMK & 12WTA"]
), "Antibody Capture .obs library_id should equal '12ABC & 12WTA."
assert "sample_tag" in prot_adata.obs.keys(), (
    "Protein modality should contain column 'sample_id'."
)
assert "sample_id" in prot_adata.obs.keys(), (
    "Protein modality should contain column 'sample_name'."
)

# check whether gene was found
assert "PDE4DIP" in data.var_names, 'Output should contain gex column "PDE4DIP".'
assert "CD279:EH12-1|PDCD1|AHS0014|pAbO" in data.var_names, (
    'Output should contain abc column "CD279:EH12-1|PDCD1|AHS0014|pAbO".'
)

print(">> Run with VDJ data", flush=True)
vdj_input = meta["resources_dir"] + "/bdrhap_vdj/sample.cellismo"
vdj_contigs = (
    meta["resources_dir"] + "/bdrhap_vdj/sample_VDJ_Unfiltered_Contigs_AIRR.tsv.gz"
)
vdj_output = "output_vdj.h5mu"

subprocess.check_output(
    [
        meta["executable"],
        "--input",
        vdj_input,
        "--input_vdj_contigs",
        vdj_contigs,
        "--output",
        vdj_output,
        "--id",
        "foo",
    ]
)
assert path.exists(vdj_output), "No output was created."

vdj_data = mu.read_h5mu(vdj_output)
assert set(vdj_data.mod.keys()) == {"rna", "vdj_t", "vdj_t_gd", "vdj_b"}, (
    f"Unexpected modalities: {list(vdj_data.mod.keys())}"
)

rna_adata = vdj_data.mod["rna"]
assert np.array_equal(
    rna_adata.var["reference_file"].unique(),
    ["BD_Rhapsody_Immune_Response_Panel_Hs.fasta"],
), "Targeted reference file should be used as reference file for targeted runs."

for modality in ["vdj_t", "vdj_t_gd", "vdj_b"]:
    vdj_adata = vdj_data.mod[modality]
    assert "airr" in vdj_adata.obsm, f"{modality} should contain .obsm['airr']."
    assert vdj_adata.n_obs > 0, f"{modality} should contain cells."
    # all contigs of all cell labels are included, not only those of putative cells
    putative = vdj_adata.obs["putative_cell"]
    assert putative.any() and not putative.all(), (
        f"{modality} should contain both putative and non-putative cells."
    )
    # putative cells correspond with the cells in the rna modality
    assert set(vdj_adata.obs_names[putative]) <= set(rna_adata.obs_names), (
        f"Putative cells of {modality} should be present in the rna modality."
    )
    assert (vdj_adata.obs["run_id"] == "foo").all()

loci = {
    "vdj_t": {"TRA", "TRB"},
    "vdj_t_gd": {"TRG", "TRD"},
    "vdj_b": {"IGH", "IGK", "IGL"},
}
for modality, expected_loci in loci.items():
    found_loci = set(ak.flatten(vdj_data.mod[modality].obsm["airr"]["locus"]))
    assert found_loci <= expected_loci, (
        f"{modality} contains unexpected loci: {found_loci - expected_loci}"
    )

print(">> Check that VDJ data without contigs file fails", flush=True)
result = subprocess.run(
    [
        meta["executable"],
        "--input",
        vdj_input,
        "--output",
        "output_vdj_fail.h5mu",
        "--id",
        "foo",
    ],
    capture_output=True,
    text=True,
    check=False,
)
assert result.returncode != 0, "Conversion should fail without --input_vdj_contigs."
assert "--input_vdj_contigs" in result.stderr, (
    "Error message should mention --input_vdj_contigs."
)

print("> Test successful", flush=True)
