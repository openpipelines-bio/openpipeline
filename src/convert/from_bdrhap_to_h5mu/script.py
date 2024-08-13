import mudata as mu

## VIASH START
par = {
    "id": "sample",
    "input": "resources_test/bdrhap_5kjrt/processed/output_raw/sample.h5mu",
    "output": "bd_rhap_to_h5mu_test.h5mu",
    "output_compression": None
}
## VIASH END

print(">> Reading input file", flush=True)
mdata = mu.read_h5mu(par["input"])

# Check if modalities are present
modalities = list(mdata.mod.keys())
assert len(modalities) > 0, "No modalities found in input data"

def process_modality_inline(adata, modality):
    adata.obs["library_id"] = " & ".join(adata.uns["Pipeline_Inputs"]["Libraries"])
    adata.obs["cell_id"] = adata.obs.index
    adata.obs["run_id"] = par["id"]
    
    adata.obs.rename(
        columns={
            "Sample_Tag": "sample_tag",
            "Sample_Name": "sample_id"},
        inplace=True)

    adata.var["gene_ids"] = adata.var.index
    adata.var["gene_name"] = adata.var.index
    
    if modality == "rna":
        adata.var["feature_type"] = "Gene Expression"
        adata.var["reference_file"] = adata.uns["Pipeline_Inputs"]["Reference_Archive"]
        
    elif modality == "prot":
        adata.var["feature_type"] = "Antibody Capture"
        adata.var["reference_file"] = " & ".join(adata.uns["Pipeline_Inputs"]["AbSeq_Reference"])
    
    # TODO: add other modalities

for key, value in mdata.mod.items():
    print(">> Processing modality:", key, flush=True)
    process_modality_inline(value, key)

print(">> Writing output file", flush=True)
mdata.write_h5mu(par["output"], compression=par["output_compression"])