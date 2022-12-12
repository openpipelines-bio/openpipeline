import mudata
import numpy as np
import scanpy as sc
import pynndescent
import numba


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "reference": "https://zenodo.org/record/6337966/files/HLCA_emb_and_metadata.h5ad",
    "targets": ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"],
    "modality": "rna",
    "reference_obsm_key": "X_integrated_scanvi",
    "query_obsm_key": "X_integrated_scanvi",
    "output": "foo.h5mu",
    "obs_output_suffix": "_pred",
    "n_neighbors": 1
}
### VIASH END


@numba.njit
def weighted_prediction(weights, ref_cats):
    """Get highest weight category."""
    N = len(weights)
    predictions = np.zeros((N,), dtype=ref_cats.dtype)
    uncertainty = np.zeros((N,))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                uncertainty[i] = max(1 - best_prob, 0)

    return predictions, uncertainty


def distances_to_affinities(distances):
    stds = np.std(distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    distances_tilda = np.exp(-np.true_divide(distances, stds))

    return distances_tilda / np.sum(distances_tilda, axis=1, keepdims=True)

def main():
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod[par["modality"]]

    adata_reference = sc.read(par["reference"], backup_url=par["reference"])

    if par["reference_obsm_key"] is None:
        X_train = adata_reference.X
    else:
        X_train = adata_reference.obsm[par["reference_obsm_key"]]

    print("X_train:")
    print(X_train)
    ref_nn_index = pynndescent.NNDescent(X_train, n_neighbors=par["n_neighbors"])
    ref_nn_index.prepare()

    if par["query_obsm_key"] is None:  # TODO: Rename to query
        query = adata.X
    else:
        query = adata.obsm[par["query_obsm_key"]]

    print("Query:")
    print(query)
    ref_neighbors, ref_distances = ref_nn_index.query(query, k=par["n_neighbors"])

    weights = distances_to_affinities(ref_distances)

    if "labels_transfer" not in adata.uns:
        adata.uns["labels_transfer"] = {}

    # for each annotation level, get prediction and uncertainty
    for target in par["targets"]:
        ref_cats = adata_reference.obs[target].cat.codes.to_numpy()[ref_neighbors]
        prediction, uncertainty = weighted_prediction(weights, ref_cats)
        prediction = np.asarray(adata_reference.obs[target].cat.categories)[prediction]
        
        predicted_label_col_name = target + par["obs_output_suffix"]
        adata.obs[predicted_label_col_name], adata.obs[target + "_uncertainty"] = prediction, uncertainty
        
        adata.uns["labels_transfer"][predicted_label_col_name] = {
            "method": "KNN_pynndescent",
            "n_neighbors": par["n_neighbors"],
            "reference": par["reference"]
        }

    mdata.mod[par['modality']] = adata
    mdata.update()
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main()
