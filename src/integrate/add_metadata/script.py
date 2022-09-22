from mudata import read_h5mu
import pandas as pd

### VIASH START
par = {
    "input": ["resources_test/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu"],
    "output": "foo.h5mu",
    "modality": ["rna"],
    "csv": ["test.csv"]
}
meta = {
    "n_proc": 10

}
### VIASH END
if par["sample_id"] and par["matrix_sample_column"]:
    raise ValueError("--sample_id and --matrix_sample_column are mutually exclusive.")
if not par["sample_id"] and not par["matrix_sample_column"]:
    raise ValueError("Must define set --sample_id or --matrix_sample_column")
metadata = pd.read_csv(par['csv'], sep=",", header=0, index_col=par["csv_sample_column"])
mdata = read_h5mu(par['input'])
for mod_name in par['modality']:
    mod_data = mdata.mod[mod_name]

    original_matrix = getattr(mod_data, par["matrix"])
    if par["sample_id"]:
        sample_metadata = metadata.loc[par['sample_id']]
        new_matrix = original_matrix.assign(**sample_metadata)
    elif par["matrix_sample_column"]:
        sample_ids = original_matrix[par['matrix_sample_column']]
        print(sample_ids.tolist())
        print(metadata)

        try:
            new_columns = metadata.loc[sample_ids.tolist()]
        except KeyError as e:
            raise KeyError(f"Not all sample IDs selected from {par['matrix']}] "
                            "(using the column selected with --matrix_sample_column) were found in "
                            "the csv file.") from e
        new_matrix = pd.concat([original_matrix.reset_index(drop=True), 
                             new_columns.reset_index(drop=True)], axis=1)\
                            .set_axis(original_matrix.index)
    setattr(mod_data, par['matrix'], new_matrix)
mdata.write_h5mu(par['output'].strip())

        

