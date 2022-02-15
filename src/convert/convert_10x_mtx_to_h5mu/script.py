import muon as mu

## VIASH START
par = {
    "input": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix",
    "output": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5mu",
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_10x_mtx(par["input"])

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
=======
    "compression": "gzip",
}
## VIASH END

mdata = mu.read_10x_mtx(par["input"])

mdata.write_h5mu(par["output"], compression=par["compression"])
>>>>>>> Add test to convert_10x_mtx_to_h5mu + correct description and example of --output + fix example input in script
=======
    "compression": "gzip",
=======
>>>>>>> Remove compression (not supported by muon) + clarify input
}
## VIASH END

mdata = mu.read_10x_mtx(par["input"])
<<<<<<< HEAD

mdata.write_h5mu(par["output"], compression=par["compression"])
>>>>>>> Add test to convert_10x_mtx_to_h5mu + correct description and example of --output + fix example input in script
=======
mdata.write_h5mu(filename=par["output"])
>>>>>>> Remove compression (not supported by muon) + clarify input
