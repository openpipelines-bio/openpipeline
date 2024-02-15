import scanpy as sc
import mudata as mu


mudata = mu.read(f'{par["input"]}')

with open(par["s_genes"]) as s_genes_file:
    s_genes = [x.strip() for x in s_genes_file]

with open(par["g2m_genes"]) as g2m_genes_file:
    g2m_genes = [x.strip() for x in g2m_genes_file]

if par["gene_pool"]:
    with open(par["gene_pool"]) as gene_pool_file:
        par["gene_pool"] = [x.strip() for x in gene_pool_file]

output = sc.tl.score_genes_cell_cycle(
    mudata.mod["rna"],
    s_genes=s_genes,
    g2m_genes=g2m_genes,
    gene_pool=par["gene_pool"],
    n_bins=par["n_bins"],
    random_state=par["random_state"],
    copy=True,
    use_raw=par["use_raw"]
)

mudata.mod["rna"] = output
mudata.write(f'{par["output"]}')
