### VIASH START
par = {
	"input": "/run/media/rcannood/Data/analyses/CS000182/output/SA00723_210416.NovaSeq1.FCB/",
  "output": "/run/media/rcannood/Data/analyses/CS000182/output/SA00723_210416.NovaSeq1.FCB/GC109678_TGCCTCTT-TATAGCCT_final.h5ad"
}
### VIASH END

import scanpy as sc
import scipy
import glob

files_in = glob.glob(par["input"] + '/*MolsPerCell_Unfiltered.csv.gz')
assert len(files_in) == 1, "Expecting only one file ending with 'MolsPerCell_Unfiltered.csv.gz' in input dir."

file_in = files_in[0]
print("Converting " + file_in + " to " + par["output"])

adata = sc.read_csv(file_in)
adata.var_names_make_unique()
adata.X = scipy.sparse.csr_matrix(adata.X)
adata.raw = adata

adata.write(par["output"], compression = "gzip")
