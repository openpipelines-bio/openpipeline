### VIASH START

par = {
}
### VIASH END

import scanpy as sc
import pandas as pd

data = sc.read_h5ad(par["input"])

sc.tl.leiden(data, resolution = float(par["resolution"]), key_added = par["clusterColumnName"])

if (par["outputFormat"] == "h5ad"):
     data.write_h5ad(par["output"], compression = "lzf")
elif (par["outputFormat"] == "csv"):
     data.obs[par["clusterColumnName"]].to_csv(par["output"])
else: 
     raise ValueError("An unrecognized output format was specified.")
