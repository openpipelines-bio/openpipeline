import mudata as mu

## VIASH START
par = {
    "id": "sample",
    "input": "resources_test/bdrhap_5kjrt/processed2/sample.h5mu",
    "output": "bd_rhap2_to_h5mu_test.h5mu",
    "output_compression": None
    } 
## VIASH END

mdata = mu.read(par["input"])
# Get modalities
modalities = list(mdata.mod.keys())
assert len(modalities) > 0, "No modalities found in input data"

# Dictionary for processed modalities
processed_modalities = {}


def process_rna_modality(adata):
    adata.obs["library_id"] = " & ".join(adata.uns["Pipeline_Inputs"]["Libraries"])
    adata.obs["cell_id"] = adata.obs.index
    adata.obs["run_id"] = par["id"]

    adata.var["gene_ids"] = adata.var.index
    adata.var["gene_name"] = adata.var.index
    adata.var["feature_type"] = "Gene Expression"
    adata.var["reference_file"] = adata.uns["Pipeline_Inputs"]["Reference_Archive"]

    return adata


def process_prot_modality(adata):
    adata.obs["library_id"] = " & ".join(adata.uns["Pipeline_Inputs"]["Libraries"])
    adata.obs["cell_id"] = adata.obs.index
    adata.obs["run_id"] = par["id"]

    adata.var["gene_ids"] = adata.var.index
    adata.var["gene_name"] = adata.var.index
    adata.var["feature_type"] = "Antibody Capture"
    adata.var["reference_file"] = " & ".join(adata.uns["Pipeline_Inputs"]["AbSeq_Reference"])

    return adata

## Processing RNA modality
if "rna" in modalities:
    rna_adata = process_rna_modality(mdata.mod["rna"])
    processed_modalities["rna"] = rna_adata


## Processing Protein modality
if "prot" in modalities:
    prot_adata = process_prot_modality(mdata.mod["prot"])
    processed_modalities["prot"] = prot_adata

##TODO: Process other modalities

output_mdata = mu.MuData(processed_modalities)
output_mdata.write_h5mu(par["output"], compression=par["output_compression"])