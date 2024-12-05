import sys
### VIASH START
par = {
    "input": "resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    "output": "test.h5mu",
    "compression": "gzip"
}
meta = {}
### VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import compress_h5mu

if __name__ == "__main__":
    compress_h5mu(par["input"], par["output"], compression=par["compression"])