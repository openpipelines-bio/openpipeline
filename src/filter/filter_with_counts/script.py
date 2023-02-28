
import mudata as mu
import numpy as np
import logging
from sys import stdout
from operator import le, ge, gt

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': 'rna',
  'output': 'output.h5mu',
  'obs_name_filter': 'filter_with_counts',
  'var_name_filter': 'filter_with_counts',
  'do_subset': True,
  'min_counts': int('200'),
  'max_counts': int('5000000'),
  'min_genes_per_cell': int('200'),
  'max_genes_per_cell': int('1500000'),
  'min_cells_per_gene': int('3'),
  'min_fraction_mito': float('0.0'),
  'max_fraction_mito': float('0.2'),
  "var_gene_names": "gene_symbol",
  "mitochondrial_gene_regex": "^[mM][tT]-"
}
meta = {
    'functionality_name': 'filter_on_counts'
}
### VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input data")
mdata = mu.read_h5mu(par["input"])

mdata.var_names_make_unique()

mod = par['modality']
logger.info("Processing modality %s.", mod)
data = mdata.mod[mod]

logger.info("\tUnfiltered data: %s", data)

logger.info("\tComputing aggregations.")
n_counts_per_cell = np.ravel(np.sum(data.X, axis=1))
n_cells_per_gene = np.sum(data.X > 0, axis=0)
n_genes_per_cell = np.sum(data.X > 0, axis=1)
genes_column = data.var[par["var_gene_names"]] if par["var_gene_names"] else data.var_names
mito_genes = genes_column.str.contains(par["mitochondrial_gene_regex"], regex=True)
pct_mito = np.ravel(np.sum(data[:, mito_genes].X, axis=1) / np.sum(data.X, axis=1))

def apply_filter_to_mask(mask, base, filter, comparator):
    new_filt = np.ravel(comparator(base, filter))
    num_removed = np.sum(np.invert(new_filt) & mask)
    mask &= new_filt
    return num_removed, mask

# Filter genes
keep_genes = np.repeat(True, data.n_vars)
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
           ("min_fraction_mito", pct_mito, ge, "\tRemoving %s cells with <%s percentage mitochondrial reads."),
           ("max_fraction_mito", pct_mito, le, "\tRemoving %s cells with >%s percentage mitochondrial reads."),
           (0, np.sum(data[:,keep_genes].X, axis=1), gt, "\tRemoving %s cells with %s counts"))

keep_cells = np.repeat(True, data.n_obs)
for filter_name_or_value, base, comparator, message in filters:
    try:
        filter = par[filter_name_or_value]
    except KeyError:
        filter = filter_name_or_value
    if filter is not None:
        num_removed, keep_cells = apply_filter_to_mask(keep_cells, base, filter, comparator)
        logger.info(message, num_removed, filter)

if par["obs_name_filter"] is not None:
    data.obs[par["obs_name_filter"]] = keep_cells
if par["var_name_filter"] is not None:
    data.var[par["var_name_filter"]] = keep_genes

if par["do_subset"]:
    mdata.mod[mod] = data[keep_cells, keep_genes]

logger.info("\tFiltered data: %s", data)
logger.info("Writing output data to %s", par["output"])
mdata.write_h5mu(par["output"], compression="gzip")

logger.info("Finished")