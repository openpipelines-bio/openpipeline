### VIASH START

par = {
}
### VIASH END

import scanpy as sc

data = sc.read_h5ad(par["input"])

sc.tl.umap(data, 
          min_dist = par["min_dist"],
          alpha = par["alpha"],
          gamma = par["gamma"],
          random_state = par["random_seed"],
          negative_sample_rate = par["negative_sample_rate"],
          init_pos = par["init_pos"])

data.write_h5ad(par["output"], compression = "gzip")
