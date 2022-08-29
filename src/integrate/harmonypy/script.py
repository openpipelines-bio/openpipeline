import mudata
from harmonypy import run_harmony


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "compression": "gzip",
    "modality": ["rna"],
    "pca_key": "X_pca",
    "output_key_suffix": "harmony",
    "theta": 0,
    "sample_key": "batch",
}
### VIASH END


def main():
    mdata = mudata.read(par["input"].strip())
    for mod_name in par['modality']:
        mod = mdata.mod[mod_name]
        pca_embedding = mod.obsm[par['pca_key']]
        metadata = mod.obs
        ho = run_harmony(pca_embedding, metadata, par['sample_key'], theta=par['theta'])
        mod.obsm[f'{par["pca_key"]}_{par["output_key_suffix"]}'] = ho.Z_corr.T
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main()