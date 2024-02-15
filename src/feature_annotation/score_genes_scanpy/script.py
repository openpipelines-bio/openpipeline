import scanpy as sc
import mudata as mu


mudata = mu.read(f'{par["input"]}')

with open(par["gene_list"]) as gene_list_file:
    gene_list = [x.strip() for x in gene_list_file]
    assert len(gene_list) > 0, "no genes detected in --gene_list"

if par["gene_pool"]:
    with open(par["gene_pool"]) as gene_pool_file:
        par["gene_pool"] = [x.strip() for x in gene_pool_file]
        assert len(par["gene_pool"]) > 0, "no genes detected in --gene_pool"

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
mudata.write(f'{par["output"]}', compression=par["output_compression"])
