import scanpy as sc
import mudata as mu

## VIASH START
par = {
    "input": "resources_test/merge_test_data/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu",
    "gene_pool": "src/annotate/score_genes_cell_cycle/test_resources/regev_lab_cell_cycle_genes.txt"
}
## VIASH END



adata = mu.read(f"{par["input"]}/rna")

output = sc.tl.score_genes(
    adata,
    par["gene_list"],
    ctrl_size = par["ctrl_size"],
    gene_pool = par["gene_pool"],
    n_bins = par["n_bins"],
    score_name = par ["score_name"],
    random_state = par["random_state"],
    copy = par["copy"],
    use_raw = par["use_raw"]
)

if output:
    mu.write(f"{par["output"]}/rna", adata)
