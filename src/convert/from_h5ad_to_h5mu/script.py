import muon as mu
import anndata

## VIASH START
par = {
    "input": ["resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"],
    "modality": ["rna"],
    "output": "output.h5mu",
    "compression": "gzip",
    "conversions_obsm": '{"counts_antibody":"prot", "counts_custom": "custom"}',
}
## VIASH END

assert len(par["input"]) == len(par["modality"]), "Number of input files should be the same length as the number of modalities"

print("Reading input files")
data = { key: anndata.read_h5ad(path) for key, path in zip(par["modality"], par["input"]) }

try:
    data.var_names_make_unique()
except:
    pass

print("Converting to muon")
muon = mu.MuData(data)

try:
    muon.var_names_make_unique()
except:
    pass

print(f"Writing {par['output']}")
muon.write_h5mu(par["output"])
