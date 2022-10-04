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
    "cpus": 10

}
### VIASH END
if par["obs_key"] and par["var_key"]:
    raise ValueError("--obs_key can not be used in conjuction with --var_key.")
if not (par["obs_key"] or par["var_key"]):
    raise ValueError("Must define either --obs_key or --var_key")
metadata = pd.read_csv(par['input_csv'], sep=",", header=0, index_col=par["csv_key"])
mdata = read_h5mu(par['input'])
mod_data = mdata.mod[par['modality']]
matrix = 'var' if par["var_key"] else 'obs'
matrix_sample_column_name = par["var_key"] if par["var_key"] else par["obs_key"]

original_matrix = getattr(mod_data, matrix)
sample_ids = original_matrix[matrix_sample_column_name]
print(sample_ids.tolist())
print(metadata)

try:
    new_columns = metadata.loc[sample_ids.tolist()]
except KeyError as e:
    raise KeyError(f"Not all sample IDs selected from {matrix} "
                    "(using the column selected with --var_key or --obs_key) were found in "
                    "the csv file.") from e
new_matrix = pd.concat([original_matrix.reset_index(drop=True), 
                        new_columns.reset_index(drop=True)], axis=1)\
                        .set_axis(original_matrix.index)    
setattr(mod_data, matrix, new_matrix)
mdata.write_h5mu(par['output'].strip())

        

