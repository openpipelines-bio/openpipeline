
import muon
import scanpy as sc
import numpy as np

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

print("Reading input data")
mdata = muon.read_h5mu(par["input"])

mdata.var_names_make_unique()

for mod in par['modality']:
    print(f"Processing modality '{mod}'")
    data = mdata.mod[mod]

    print(f"  Unfiltered data: {data}")

    print("  Computing aggregations")
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
        print(f"  Removing {num_removed} genes with non-zero values in <{par['min_cells_per_gene']} cells.")
        keep_genes &= new_filt

    if par["min_genes_per_cell"] is not None:
        new_filt = np.ravel(n_genes_per_cell >= par['min_genes_per_cell'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with non-zero values in <{par['min_genes_per_cell']} genes.")
        keep_cells &= new_filt

    if par["max_genes_per_cell"] is not None:
        new_filt = np.ravel(n_genes_per_cell <= par['max_genes_per_cell'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with non-zero values in >{par['max_genes_per_cell']} genes.")
        keep_cells &= new_filt

    if par["min_counts"] is not None:
        new_filt = np.ravel(n_counts_per_cell >= par['min_counts'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with <{par['min_counts']} total counts.")
        keep_cells &= new_filt

    if par["max_counts"] is not None:
        new_filt = np.ravel(n_counts_per_cell <= par['max_counts'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with >{par['max_counts']} total counts.")
        keep_cells &= new_filt

    if par["min_fraction_mito"] is not None:
        new_filt = np.ravel(pct_mito >= par['min_fraction_mito'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with <{par['min_fraction_mito']} percentage mitochondrial realds.")
        keep_cells &= new_filt

    if par["max_fraction_mito"] is not None:
        new_filt = np.ravel(pct_mito <= par['max_fraction_mito'])
        num_removed = np.sum(np.invert(new_filt) & keep_cells)
        print(f"  Removing {num_removed} cells with >{par['max_fraction_mito']} percentage mitochondrial realds.")
        keep_cells &= new_filt
    
    # remove cells with zero counts
    new_filt = np.ravel(np.sum(data[:,keep_genes].X, axis=1)) > 0
    num_removed = np.sum(np.invert(new_filt) & keep_cells)
    print(f"  Removing {num_removed} cells with zero counts")
    keep_cells &= new_filt

    if par["obs_name_filter"] is not None:
        data.obs[par["obs_name_filter"]] = keep_cells
    if par["var_name_filter"] is not None:
        data.var[par["var_name_filter"]] = keep_genes

    if par["do_subset"]:
        mdata.mod[mod] = data[keep_cells, keep_genes]
    
    print(f"  Filtered data: {data}")

# can we assume execution_log exists?
if mdata.uns is None or "execution_log" not in mdata.uns:
    mdata.uns["execution_log"] = []
# store new entry
new_entry = {"component": meta["functionality_name"], "params": par}
mdata.uns["execution_log"] = str(mdata.uns["execution_log"] + [new_entry])


print("Writing output data")
mdata.write(par["output"])
