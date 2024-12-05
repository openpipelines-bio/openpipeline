import mudata

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "bar.h5mu",
    "number_of_observations": 100,
}
### VIASH END

if __name__ == "__main__":
    # read data
    data = mudata.read(par["input"])

    # subset data
    if par["modality"]:
        data.mod[par["modality"]] = data.mod[par["modality"]][
            : par["number_of_observations"]
        ]
    else:
        data = data[: par["number_of_observations"]]

    # write data
    data.write_h5mu(par["output"], compression=par["output_compression"])
