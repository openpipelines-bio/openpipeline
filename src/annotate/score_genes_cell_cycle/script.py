import scanpy as sc
import mudata as mu

## VIASH START
par = {
    "input": 
}
## VIASH END

adata = mu.read(f"{par["input"]}/rna")

output = sc.tl.score_genes_cell_cycle(
    adata,
    s_genes = par["s_genes"],
    g2m_genes = par["g2m_genes"],
    gene_pool = par["gene_pool"],
    n_bins = par["n_bins"],
    score_name = par ["score_name"],
    random_state = par["random_state"],
    copy = par["copy"],
    use_raw = par["use_raw"]
)

if output:
    mu.write(f"{par["output"]}/rna", adata)