import scanpy as sc
import mudata as mu


mudata = mu.read(f'{par["input"]}')

par["use_raw"] = par.setdefault("use_raw", None)
par["gene_pool"] = par.setdefault("gene_pool", None)

gene_list = [x.strip() for x in open(par["gene_list"])]

if par["gene_pool"]:
    par["gene_pool"] = [x.strip() for x in open(par["gene_pool"])]

output = sc.tl.score_genes(
    mudata.mod["rna"],
    gene_list,
    ctrl_size=par["ctrl_size"],
    gene_pool=par["gene_pool"],
    n_bins=par["n_bins"],
    score_name=par["score_name"],
    random_state=par["random_state"],
    copy=True,
    use_raw=par["use_raw"]
)

mudata.mod["rna"] = output
mudata.write(f'{par["output"]}')
