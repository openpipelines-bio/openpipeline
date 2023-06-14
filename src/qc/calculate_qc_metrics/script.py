from mudata import read_h5mu
from scipy.sparse import issparse, isspmatrix_coo, csr_matrix
from sklearn.utils.sparsefuncs import mean_variance_axis
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo.h5mu",
    "modality": "rna",
    "layer": None,
    "top_n_vars": [10, 20, 50],
    "var_qc_metrics": None,
}
## VIASH END

def main():
    input_data = read_h5mu(par["input"])
    modality_data = input_data.mod[par["modality"]]
    var = modality_data.var
    layer = modality_data.X if not par['layer'] else modality_data.layers[par['layer']]
    if not issparse(layer):
        raise NotImplementedError("Expected layer to be in sparse format.")
    if isspmatrix_coo(layer):
        layer = csr_matrix(layer)
    layer.eliminate_zeros()
    
    # var statistics
    num_nonzero_obs = layer.getnnz(axis=0)
    obs_mean, _  = mean_variance_axis(layer, axis=0)
    pct_dropout = (1 - num_nonzero_obs / layer.shape[0]) * 100
    total_counts_obs = np.ravel(layer.sum(axis=0))

    # obs statistics
    num_nonzero_vars = layer.getnnz(axis=1)
    total_counts_var = np.ravel(layer.sum(axis=1))

    top_metrics = {}
    if par["top_n_vars"]:
        par["top_n_vars"] = sorted(par["top_n_vars"])
        distributions = get_top_from_csr_matrix(layer, par["top_n_vars"])
        top_metrics = {distribution_size: distribution * 100
                       for distribution_size, distribution 
                       in zip(par["top_n_vars"], distributions.T)}
    
    total_expr_qc = {}
    pct_expr_qc = {}
    if par["var_qc_metrics"]:
        for qc_metric in par["var_qc_metrics"]:
            if not qc_metric in var:
                raise ValueError(f"Value for --var_qc_metrics, {qc_metric} "
                                 f"not found in .var for modality {par['modality']}")
            qc_column = var[qc_metric].values
            if not set(np.unique(qc_column)) == set((True, False)):
                raise ValueError(f"Column {qc_metric} in .var for modality {par['modality']} "
                                 f"must only contain boolean values")
            
            total_expr_qc[qc_metric] = np.ravel(layer[:, qc_column].sum(axis=1))
            pct_expr_qc[qc_metric] =  total_expr_qc[qc_metric] / total_counts_var * 100
    
    # Write all of the calculated statistics
    modality_data.var = modality_data.var.assign(
        **{"pct_dropout": pct_dropout,
           "num_nonzero_obs": num_nonzero_obs,
           "obs_mean": obs_mean,
           "total_counts": total_counts_obs})
    
    modality_data.obs = modality_data.obs.assign(
        **({"num_nonzero_vars": num_nonzero_vars,
            "total_counts": total_counts_var} | \
           {f"pct_{qc_metric}": col for qc_metric, col in pct_expr_qc.items()} | \
           {f"total_counts_{qc_metrix}": col for qc_metrix, col in total_expr_qc.items()}) | \
           {f"pct_of_counts_in_top_{n_top}_vars": col for n_top, col in top_metrics.items()})

    input_data.write(par["output"], compression=par["output_compression"])
            
def get_top_from_csr_matrix(matrix, top_n_genes):
    # csr matrices stores a 3D matrix in a format such that data for individual cells
    # are stored in 1 array. Another array (indptr) here stores the ranges of indices
    # to select from the data-array (.e.g. data[indptr[0]:indptr[1]] for row 0) for each row.
    # Another array 'indices' maps each element of data to a column 
    # (data and indices arrays have the same length)
    top_n_genes = np.array(top_n_genes).astype(np.int64)
    assert np.all(top_n_genes[:-1] <= top_n_genes[1:]), "top_n_genes must be sorted"
    row_indices, data = matrix.indptr, matrix.data
    number_of_rows, max_genes_to_parse = row_indices.size-1, top_n_genes[-1]
    top_data = np.zeros((number_of_rows, max_genes_to_parse), 
                        dtype=data.dtype)
    # Loop over each row to create a dense matrix without the 0 counts,
    # but not for the whole matrix, only store the genes up until
    # the largest number of top n genes.
    for row_number in range(number_of_rows):
        row_start_index, row_end_index = row_indices[row_number], row_indices[row_number+1]
        row_data = data[row_start_index:row_end_index] # all non-zero counts for an row
        try:
            # There are less genes with counts in the row than the
            # maximum number of genes we would like to select
            # all these genes are in the top genes, just store them
            top_data[row_number, :row_end_index-row_start_index] = row_data
        except ValueError:
            # Store the counts for the top genes
            top_data[row_number, :] = np.partition(row_data, -max_genes_to_parse)[-max_genes_to_parse:]

    # Partition works from smallest to largest, but we want largest
    # so do smallest to largest first (but with reversed indices)
    top_data = np.partition(top_data, max_genes_to_parse - top_n_genes)
    # And then switch the order around
    top_data = np.flip(top_data, axis=1)

    cumulative = top_data.cumsum(axis=1, dtype=np.float64)[:,top_n_genes-1]
    return cumulative / np.array(matrix.sum(axis=1))

if __name__ == "__main__":
    main()