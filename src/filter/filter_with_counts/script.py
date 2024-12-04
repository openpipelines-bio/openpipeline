
import mudata as mu
import numpy as np
import sys
from operator import le, ge, gt

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': 'rna',
  'output': 'output.h5mu',
  'obs_name_filter': 'filter_with_counts',
  'var_name_filter': 'filter_with_counts',
  'do_subset': True,
  'min_counts': 200,
  'max_counts': 5000000,
  'min_genes_per_cell': 200,
  'max_genes_per_cell': 1500000,
  'min_cells_per_gene': 3,
  'layer': None
}
meta = {
    'name': 'filter_on_counts',
    'resources_dir': '.'
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

logger.info("Reading input data")
mdata = mu.read_h5mu(par["input"])

mdata.var_names_make_unique()

mod = par['modality']
logger.info("Processing modality %s.", mod)
modality_data = mdata.mod[mod]
logger.info("\tUnfiltered data: %s", modality_data)

logger.info("Selecting input layer %s", "X" if par["layer"] else par["layer"])
input_layer = modality_data.X if not par["layer"] else modality_data.layers[par["layer"]]

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
    num_removed, keep_genes = apply_filter_to_mask(keep_genes,
                                                   n_cells_per_gene,
                                                   par['min_cells_per_gene'],
                                                   ge)
    logger.info("\tRemoving %s genes with non-zero values in <%s cells.",
                num_removed, par['min_cells_per_gene'])

# Filter cells
filters = (("min_genes_per_cell", n_genes_per_cell, ge, "\tRemoving %s cells with non-zero values in <%s genes."),
           ("max_genes_per_cell", n_genes_per_cell, le, "\tRemoving %s cells with non-zero values in >%s genes."),
           ("min_counts", n_counts_per_cell, ge, "\tRemoving %s cells with <%s total counts."),
           ("max_counts", n_counts_per_cell, le, "\tRemoving %s cells with >%s total counts."),
           (0, np.sum(input_layer[:,keep_genes], axis=1), gt, "\tRemoving %s cells with %s counts"))

keep_cells = np.repeat(True, modality_data.n_obs)
for filter_name_or_value, base, comparator, message in filters:
    try:
        filter = par[filter_name_or_value]
    except KeyError:
        filter = filter_name_or_value
    if filter is not None:
        num_removed, keep_cells = apply_filter_to_mask(keep_cells, base, filter, comparator)
        logger.info(message, num_removed, filter)

if par["obs_name_filter"] is not None:
    modality_data.obs[par["obs_name_filter"]] = keep_cells
if par["var_name_filter"] is not None:
    modality_data.var[par["var_name_filter"]] = keep_genes

if par["do_subset"]:
    mdata.mod[mod] = modality_data[keep_cells, keep_genes]

logger.info("\tFiltered data: %s", modality_data)
logger.info("Writing output data to %s", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])

logger.info("Finished")