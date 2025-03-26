from umap import UMAP
import mudata as mu
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "obsm_output": "X_densmap",
    "lambda": 2.0,
    "fraction": 0.3,
    "var_shift": 0.1,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading %s, modality %s", par["input"], par["modality"])
try:
    data = mu.read_h5ad(par["input"], mod=par["modality"])
except KeyError as e:
    raise ValueError(
        f"Modality '{par['modality']}' not found in the input data."
    ) from e

logger.info("Computing densMAP for modality '%s'", par["modality"])

neigh_key = par["uns_neighbors"]

if neigh_key not in data.uns:
    raise ValueError(
        f"'{neigh_key}' was not found in .mod['{par['modality']}'].uns. Set the correct key or run 'find_neighbors' first."
    )

temp_uns = {neigh_key: data.uns[neigh_key]}

if "use_rep" not in temp_uns[neigh_key]["params"]:
    raise ValueError(
        f"'use_rep' was not found in .mod['{par['modality']}'].uns['{neigh_key}'].params. Set the correct key or run PCA first."
    )


X_densmap = UMAP(
    min_dist=par["min_dist"],
    spread=par["spread"],
    n_components=par["num_components"],
    n_epochs=par["max_iter"],
    learning_rate=par["alpha"],
    repulsion_strength=par["gamma"],
    negative_sample_rate=par["negative_sample_rate"],
    init=par["init_pos"],
    metric=data.uns["neighbors"].get("metric", "euclidean"),
    metric_kwds=data.uns["neighbors"].get("metric_kwds", {}),
    densmap=True,
    dens_lambda=par["lambda"],
    dens_frac=par["fraction"],
    dens_var_shift=par["var_shift"],
).fit_transform(data.obsm[par["obsm_pca"]])

logger.info(
    f"Writing densMAP embeddings to .mod[{par['modality']}].obsm[{par['obsm_output']}]"
)
data.obsm[par["obsm_output"]] = X_densmap

logger.info(f"Writing densMAP metadata to .mod[{par['modality']}].uns['densmap']")
data.uns["densmap"] = {
    "params": {
        "min_dist": par["min_dist"],
        "spread": par["spread"],
        "n_components": par["num_components"],
        "n_epochs": par["max_iter"],
        "learning_rate": par["alpha"],
        "repulsion_strength": par["gamma"],
        "negative_sample_rate": par["negative_sample_rate"],
        "init": par["init_pos"],
        "metric": data.uns["neighbors"].get("metric", "euclidean"),
        "metric_kwds": data.uns["neighbors"].get("metric_kwds", {}),
        "dens_lambda": par["lambda"],
        "dens_frac": par["fraction"],
        "dens_var_shift": par["var_shift"],
    }
}

logger.info("Writing to %s.", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)

logger.info("Finished")
