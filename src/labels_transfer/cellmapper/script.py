import sys

## VIASH START
par = {
    "input": "query.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    "reference": "reference.h5mu",
    "reference_obsm_features": "X_pca",
    "reference_obs_targets": ["cell_type"],
    "output": "output.h5mu",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_compression": "gzip",
    "n_neighbors": 15,
    "kernel_method": None,
}
meta = {"resources_dir": "src/labels_transfer/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import check_arguments


def main(par):
    check_arguments(par)


if __name__ == "__main__":
    sys.exit(main(par))
