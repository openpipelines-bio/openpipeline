import scanpy as sc
import muon as mu

## VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu',
  'modality': ['rna'],
  'output': 'output.h5mu',
  'output_key': 'umap',
  'min_dist': 0.5,
  'spread': 1.0,
  'num_components': 2,
  'max_iter': None,
  'alpha': 1.0,
  'gamma': 1.0,
  'negative_sample_rate': 5,
  'init_pos': 'spectral'
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    print(f"Computing UMAP for modality '{mod}'")
    data = mdata.mod[mod]

    sc.tl.umap(
        data,
        min_dist=par["min_dist"],
        spread=par["spread"],
        n_components=par["num_components"],
        maxiter=par["max_iter"],
        alpha=par["alpha"],
        gamma=par["gamma"],
        negative_sample_rate=par["negative_sample_rate"],
        init_pos=par["init_pos"],
    )
    # note: should be able to set the neighbors key

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
