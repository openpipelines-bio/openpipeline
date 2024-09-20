import mudata as mu
import anndata as ad

## VIASH START
par = {
  "input_dataset": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  "input_output": "dr_method_output.h5ad",
  "input_modality": "rna",
  "output_obsm_key": "X_dr",
  "output": "output.h5mu"
}
## VIASH END

print("Reading h5mu file", flush=True)
mdata = mu.read_h5mu(par["input_dataset"])
adata = ad.read_h5ad(par["input_output"])

mdata.mod[par["input_modality"]].obsm[par["output_obsm_key"]] = adata.obsm["X_emb"]

print("Writing h5ad file", flush=True)
mdata.write_h5mu(par["output"])
