import sys

import h5py
import scanpy as sc

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

logger.info("Reading '%s'", par["input"])
adata = sc.read_10x_h5(par["input"], gex_only=False)

logger.info("Filtering out cells with library size < %d", par["min_library_size"])
cell_mask, _ = sc.pp.filter_cells(
    adata, min_counts=par["min_library_size"], inplace=False
)

logger.info("Filtering out features with num cells < %d", par["min_cells_per_gene"])
gene_mask, _ = sc.pp.filter_genes(
    adata[cell_mask], min_cells=par["min_cells_per_gene"], inplace=False
)
keep_types = par["keep_feature_types"] or []
feature_mask = gene_mask | adata.var["feature_types"].isin(keep_types).to_numpy()

# 10x h5 files store a features x barcodes matrix in CSC format
mat = adata[cell_mask, feature_mask].X.T.tocsc()
mat.sort_indices()
logger.info("Kept %d cells and %d features", mat.shape[1], mat.shape[0])

# scanpy cannot write 10x h5 files, so write the output with h5py,
# using the same dtypes and file attributes as the input.
logger.info("Writing '%s'", par["output"])
with (
    h5py.File(par["input"], "r") as h5_in,
    h5py.File(par["output"], "w") as h5_out,
):

    def write(path, value):
        dtype = h5_in[path].dtype
        h5_out.create_dataset(path, data=value, dtype=dtype, compression="gzip")

    h5_out.attrs.update(h5_in.attrs)
    write("matrix/data", mat.data)
    write("matrix/indices", mat.indices)
    write("matrix/indptr", mat.indptr)
    write("matrix/shape", mat.shape)
    write("matrix/barcodes", h5_in["matrix/barcodes"][()][cell_mask])
    for key, dataset in h5_in["matrix/features"].items():
        value = dataset[()]
        # _all_tag_keys lists the names of extra feature fields, it is not per feature
        if key != "_all_tag_keys":
            value = value[feature_mask]
        write(f"matrix/features/{key}", value)
