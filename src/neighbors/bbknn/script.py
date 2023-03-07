from mudata import read_h5mu
import bbknn

### VIASH START
par = {
    'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu',
    'modality': 'rna',
    'obs_batch': 'sample_id',
    'obsm_input': 'X_pca',
    'n_neighbors_within_batch': 3,
    'n_trim': None,
    'n_pcs': 50
}
### VIASH END

# todo: make use of --uns_output, --obsp_connectivities and --obsp_distances
# to configure output field

h5mu_data = read_h5mu(par["input"])
modality_name = par["modality"]
modality = h5mu_data.mod[modality_name]
bbknn.bbknn(
    modality,
    use_rep=par["obsm_input"],
    batch_key = par["obs_batch"],
    neighbors_within_batch=par["n_neighbors_within_batch"],
    n_pcs=par["n_pcs"],
    trim=par["n_trim"]
)

h5mu_data.write_h5mu(par["output"], compression=par["output_compression"])