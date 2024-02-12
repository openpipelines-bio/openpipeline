import scanpy as sc
import mudata as mu


mudata = mu.read(f'{par["input"]}')

par["use_raw"] = par.setdefault("use_raw", None)
par["gene_pool"] = par.setdefault("gene_pool", None)

s_genes = [x.strip() for x in open(par["s_genes"])]
g2m_genes = [x.strip() for x in open(par["g2m_genes"])]

if par["gene_pool"]:
    par["gene_pool"] = [x.strip() for x in open(par["gene_pool"])]

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
