import scvelo as scv

## VIASH START
par = {
  'input': 'output.loom',
  'output': 'output.h5ad'
}
## VIASH END


adata = scv.read(par['input'])

adata_orig = adata.copy()

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

adata.write_h5ad(par['output'])

# TODO: convert to muon

# scv.pl.velocity_embedding_stream(adata, basis='umap')
