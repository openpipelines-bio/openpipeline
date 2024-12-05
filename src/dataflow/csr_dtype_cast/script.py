import mudata as mu
from scipy.sparse import issparse
import numpy as np

## VIASH START
par = {
    "input": "test.h5mu",
    "dtype": "int32",
    "output_compression": None,
    "modality": "rna",
    "layer": None,
    "output": "test_casted.h5mu"
}

import anndata as ad
import pandas as pd
from scipy.sparse import random
rng = np.random.default_rng(seed=1)

random_counts = random(
    50000, 100,
    density=0.8,
    format='csr',
    dtype=np.uint32,
    random_state=rng
    )
bad_dtype = random_counts.astype(np.float32)
bad_dtype.indptr = bad_dtype.indptr.astype(np.uint32)
bad_dtype.indices = bad_dtype.indices.astype(np.int64)
del random_counts 
mod1 = ad.AnnData(
    X=bad_dtype,
    obs=pd.DataFrame(index=pd.RangeIndex(50000)),
    var=pd.DataFrame(index=pd.RangeIndex(100))
    )

mdata = mu.MuData({"mod1": mod1})
mdata.write(par["input"])
par["modality"] = "mod1"
## VIASH END

# Read in data
mdata = mu.read_h5mu(par["input"])
adata = mdata.mod[par["modality"]]
layer = adata.X if not par['layer'] else adata.layers[par['layer']]
if not issparse(layer):
    raise NotImplementedError("Expected layer to be in sparse format.")

# cast dtypes of indptr and indices of csr to scipy-compatible dtypes
if not layer.indices.dtype == par["dtype"]:
    layer.indices = layer.indices.astype(par["dtype"])
if not layer.indptr.dtype == par["dtype"]:
    layer.indptr = layer.indptr.astype(par["dtype"])

# write data with casted dtype layers back to file
if not par["layer"]:
    adata.X = layer
else:
    adata.layers[par["layer"]] = layer

mdata.write_h5mu(par["output"], compression=par["output_compression"])
