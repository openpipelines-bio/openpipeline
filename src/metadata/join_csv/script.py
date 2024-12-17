import sys
import pandas as pd
from mudata import read_h5mu

### VIASH START
par = {
    "input": "work/f5/5f6365898ca5a42a360301a0c9e200/TSP15_Eye_ScleraEtc_10X_2_1.add_id.output.h5mu",
    "input_csv": "work/f5/5f6365898ca5a42a360301a0c9e200/sample_info.csv",
    "output": "foo.h5mu",
    "modality": "rna",
    "csv_key": "id",
    "obs_key": "sample_id",
    "var_key": None,
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

if par["obs_key"] and par["var_key"]:
    raise ValueError("--obs_key can not be used in conjuction with --var_key.")
if not (par["obs_key"] or par["var_key"]):
    raise ValueError("Must define either --obs_key or --var_key")

logger.info("Read metadata csv from file")
metadata = pd.read_csv(par["input_csv"], sep=",", header=0, index_col=par["csv_key"])
metadata.fillna("", inplace=True)

logger.info("Read mudata from file")
mdata = read_h5mu(par["input"])
mod_data = mdata.mod[par["modality"]]

logger.info("Joining csv to mudata")
matrix = "var" if par["var_key"] else "obs"
matrix_sample_column_name = par["var_key"] if par["var_key"] else par["obs_key"]
original_matrix = getattr(mod_data, matrix)
sample_ids = original_matrix[matrix_sample_column_name]

try:
    new_columns = metadata.loc[sample_ids.tolist()]
except KeyError as e:
    raise KeyError(
        f"Not all sample IDs selected from {matrix} "
        "(using the column selected with --var_key or --obs_key) were found in "
        "the csv file."
    ) from e
new_matrix = pd.concat(
    [original_matrix.reset_index(drop=True), new_columns.reset_index(drop=True)], axis=1
).set_axis(original_matrix.index)
setattr(mod_data, matrix, new_matrix)

logger.info("Write output to mudata file")
mdata.write_h5mu(par["output"], compression=par["output_compression"])
