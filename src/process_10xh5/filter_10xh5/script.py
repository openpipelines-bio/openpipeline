import sys

import h5py
import numpy as np
import scipy.sparse as sp

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
    "output": "output.h5",
    "min_library_size": 1000,
    "min_cells_per_gene": 300,
    "keep_feature_types": ["Antibody Capture"],
    "verbose": True,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
if not par["verbose"]:
    logger.setLevel("WARNING")

with (
    h5py.File(par["input"], "r") as h5_in,
    h5py.File(par["output"], "w") as h5_out,
):
    # 10x h5 files store a features x barcodes matrix in CSC format.
    # Read only the matrix itself, keeping the original integer dtype.
    logger.info("Reading '%s'", par["input"])
    matrix = h5_in["matrix"]
    mat = sp.csc_matrix(
        (matrix["data"][()], matrix["indices"][()], matrix["indptr"][()]),
        shape=tuple(matrix["shape"][()]),
    )

    logger.info("Filtering out cells with library size < %d", par["min_library_size"])
    library_size = np.asarray(mat.sum(axis=0, dtype=np.int64)).ravel()
    cell_mask = library_size >= par["min_library_size"]
    mat = mat[:, cell_mask]

    logger.info("Filtering out features with num cells < %d", par["min_cells_per_gene"])
    cells_per_feature = np.bincount(mat.indices[mat.data > 0], minlength=mat.shape[0])
    gene_mask = cells_per_feature >= par["min_cells_per_gene"]
    keep_types = par["keep_feature_types"] or []
    feature_types = matrix["features/feature_type"][()].astype(str)
    feature_mask = gene_mask | np.isin(feature_types, keep_types)

    mat = mat[feature_mask, :].tocsc()
    mat.sort_indices()
    logger.info("Kept %d cells and %d features", mat.shape[1], mat.shape[0])

    # Write the output using the same dtypes and file attributes as the input.
    logger.info("Writing '%s'", par["output"])

    def write(path, value):
        dtype = h5_in[path].dtype
        h5_out.create_dataset(path, data=value, dtype=dtype, compression="gzip")

    h5_out.attrs.update(h5_in.attrs)
    write("matrix/data", mat.data)
    write("matrix/indices", mat.indices)
    write("matrix/indptr", mat.indptr)
    write("matrix/shape", mat.shape)
    write("matrix/barcodes", h5_in["matrix/barcodes"][()][cell_mask])
    # Position of each kept feature in the output, used to update feature indices
    new_feature_index = np.cumsum(feature_mask) - 1
    for key, item in h5_in["matrix/features"].items():
        path = f"matrix/features/{key}"
        if key == "target_sets":
            # Probe-based matrices (Flex, Targeted GEX, Visium FFPE) have one
            # dataset per panel, listing the indices of the targeted features
            # (stored as strings by Cell Ranger).
            for panel, dataset in item.items():
                indices = dataset[()].astype(np.int64)
                indices = indices[feature_mask[indices]]
                new_indices = new_feature_index[indices].astype(dataset.dtype)
                write(f"{path}/{panel}", new_indices)
        elif isinstance(item, h5py.Group):
            raise ValueError(f"Unsupported group '{path}' in '{par['input']}'")
        elif key == "_all_tag_keys":
            # Lists the names of extra feature fields, it is not per feature
            write(path, item[()])
        else:
            write(path, item[()][feature_mask])
