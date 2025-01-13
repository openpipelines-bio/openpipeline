import mudata
from harmonypy import run_harmony


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "compression": "gzip",
    "modality": "rna",
    "obsm_input": "X_pca",
    "obsm_output": "X_pca_harmonypy",
    "theta": 2,
    "obs_covariates": ["sample_id"],
}
### VIASH END


def main():
    mdata = mudata.read(par["input"].strip())
    mod_name = par["modality"]
    mod = mdata.mod[mod_name]
    pca_embedding = mod.obsm[par["obsm_input"]]
    metadata = mod.obs
    ho = run_harmony(pca_embedding, metadata, par["obs_covariates"], theta=par["theta"])
    mod.obsm[par["obsm_output"]] = ho.Z_corr.T
    mdata.write_h5mu(par["output"].strip(), compression=par["output_compression"])


if __name__ == "__main__":
    main()
