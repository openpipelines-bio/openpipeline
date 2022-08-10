
import muon
import numpy as np
import logging
from sys import stdout

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': [ 'rna' ],
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
  'max_fraction_mito': float('0.2')
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
mdata = muon.read_h5mu(par["input"])

mdata.var_names_make_unique()

for mod in par['modality']:
    logger.info("Processing modality %s.", mod)
    data = mdata.mod[mod]

    logger.info("\tUnfiltered data: %s", data)

    logger.info("\tComputing aggregations.")
    n_counts_per_cell = np.ravel(np.sum(data.X, axis=1))
    n_cells_per_gene = np.sum(data.X > 0, axis=0)
    n_genes_per_cell = np.sum(data.X > 0, axis=1)
    keep_cells = np.repeat(True, data.n_obs)
    keep_genes = np.repeat(True, data.n_vars)

    mito_genes = data.var_names.str.contains("^[mM][tT]-")
    pct_mito = np.ravel(np.sum(data[:, mito_genes].X, axis=1) / np.sum(data.X, axis=1))

    if par["min_cells_per_gene"] is not None:
        new_filt = np.ravel(n_cells_per_gene >= par['min_cells_per_gene'])
        num_removed = np.sum(np.invert(new_filt) & keep_genes)
        logger.info("\tRemoving %s genes with non-zero values in <%s cells.", num_removed, par['min_cells_per_gene'])
        keep_genes &= new_filt

    if par["min_genes_per_cell"] is not None:
        new_filt = np.ravel(n_genes_per_cell >= par['min_genes_per_cell'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with non-zero values in <%s genes.", num_removed, par['min_genes_per_cell'])
        keep_cells &= new_filt

    if par["max_genes_per_cell"] is not None:
        new_filt = np.ravel(n_genes_per_cell <= par['max_genes_per_cell'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with non-zero values in >%s genes.", num_removed, par['max_genes_per_cell'])
        keep_cells &= new_filt

    if par["min_counts"] is not None:
        new_filt = np.ravel(n_counts_per_cell >= par['min_counts'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with <%s total counts.", num_removed, par['min_counts'])
        keep_cells &= new_filt

    if par["max_counts"] is not None:
        new_filt = np.ravel(n_counts_per_cell <= par['max_counts'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with >%s total counts.", num_removed, par['max_counts'])
        keep_cells &= new_filt

    if par["min_fraction_mito"] is not None:
        new_filt = np.ravel(pct_mito >= par['min_fraction_mito'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with <%s percentage mitochondrial reads.", num_removed, par['min_fraction_mito'])
        keep_cells &= new_filt

    if par["max_fraction_mito"] is not None:
        new_filt = np.ravel(pct_mito <= par['max_fraction_mito'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        logger.info("\tRemoving %s cells with >%s percentage mitochondrial reads.", num_removed, par['max_fraction_mito'])
        keep_cells &= new_filt
    
    # remove cells with zero counts
    new_filt = np.ravel(np.sum(data[:,keep_genes].X, axis=1)) > 0
    num_removed = np.sum(np.invert(new_filt) & keep_cells)
    logger.info("\tRemoving %s cells with zero counts", num_removed)
    keep_cells &= new_filt

    if par["obs_name_filter"] is not None:
        data.obs[par["obs_name_filter"]] = keep_cells
    if par["var_name_filter"] is not None:
        data.var[par["var_name_filter"]] = keep_genes

    if par["do_subset"]:
        mdata.mod[mod] = data[keep_cells, keep_genes]
    
    logger.info("\tFiltered data: %s", data)

logger.info("Writing output data to %s", par["output"])
mdata.write(par["output"])

logger.info("Finished")