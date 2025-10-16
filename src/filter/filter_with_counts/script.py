import mudata as mu
import numpy as np
import sys
from operator import le, ge, gt
import scipy

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "obs_name_filter": "filter_with_counts",
    "var_name_filter": "filter_with_counts",
    "do_subset": True,
    "min_counts": 200,
    "max_counts": 5000000,
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 1500000,
    "min_cells_per_gene": 3,
    "layer": None,
}
meta = {"name": "filter_on_counts", "resources_dir": "."}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading input data from %s, modality %s", par["input"], par["modality"])
modality_data = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("\tUnfiltered data: %s", modality_data)

layer_name = "X" if not par["layer"] else par["layer"]
logger.info("Selecting input layer %s", layer_name)
input_layer = (
    modality_data.X if not par["layer"] else modality_data.layers[par["layer"]]
)
if scipy.sparse.issparse(input_layer):
    # Below is a check that scipy does not do with check_format; see https://github.com/scipy/scipy/issues/23784
    # also, check_format might adjust the matrix attributes in place; so do this check before check_format

    assert input_layer.nnz == input_layer.data.size, (
        f"The provided sparse matrix (i.e. {layer_name} from modality {par['modality']}) is malformatted. "
        f"The number of elements that the matrix should hold (i.e. .nnz={input_layer.nnz}) does not correspond with "
        f"the number of elements actually being stored (i.e. .data.size={input_layer.data.size})"
    )
    input_layer.check_format(full_check=True)


logger.info("\tComputing aggregations.")
n_counts_per_cell = np.ravel(np.sum(input_layer, axis=1))
n_cells_per_gene = np.sum(input_layer > 0, axis=0)
n_genes_per_cell = np.sum(input_layer > 0, axis=1)


def apply_filter_to_mask(mask, base, filter, comparator):
    new_filt = np.ravel(comparator(base, filter))
    num_removed = np.sum(np.invert(new_filt) & mask)
    mask &= new_filt
    return num_removed, mask


# Filter genes
keep_genes = np.repeat(True, modality_data.n_vars)
if par["min_cells_per_gene"] is not None:
    num_removed, keep_genes = apply_filter_to_mask(
        keep_genes, n_cells_per_gene, par["min_cells_per_gene"], ge
    )
    logger.info(
        "\tRemoving %s genes with non-zero values in <%s cells.",
        num_removed,
        par["min_cells_per_gene"],
    )

# Filter cells
filters = (
    (
        "min_genes_per_cell",
        n_genes_per_cell,
        ge,
        "\tRemoving %s cells with non-zero values in <%s genes.",
    ),
    (
        "max_genes_per_cell",
        n_genes_per_cell,
        le,
        "\tRemoving %s cells with non-zero values in >%s genes.",
    ),
    ("min_counts", n_counts_per_cell, ge, "\tRemoving %s cells with <%s total counts."),
    ("max_counts", n_counts_per_cell, le, "\tRemoving %s cells with >%s total counts."),
    (
        0,
        np.sum(input_layer[:, keep_genes], axis=1),
        gt,
        "\tRemoving %s cells with %s counts",
    ),
)

keep_cells = np.repeat(True, modality_data.n_obs)
for filter_name_or_value, base, comparator, message in filters:
    try:
        filter = par[filter_name_or_value]
    except KeyError:
        filter = filter_name_or_value
    if filter is not None:
        num_removed, keep_cells = apply_filter_to_mask(
            keep_cells, base, filter, comparator
        )
        logger.info(message, num_removed, filter)

if par["obs_name_filter"] is not None:
    modality_data.obs[par["obs_name_filter"]] = keep_cells
if par["var_name_filter"] is not None:
    modality_data.var[par["var_name_filter"]] = keep_genes

if par["do_subset"]:
    modality_data = modality_data[keep_cells, keep_genes]

logger.info("\tFiltered data: %s", modality_data)
logger.info(
    "Writing output data to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"],
    par["input"],
    par["modality"],
    modality_data,
    par["output_compression"],
)


logger.info("Finished")
