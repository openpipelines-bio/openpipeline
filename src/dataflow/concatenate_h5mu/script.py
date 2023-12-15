from __future__ import annotations
import sys
import anndata
import mudata as mu
import pandas as pd
import numpy as np
from collections.abc import Iterable
from multiprocessing import Pool
from pathlib import Path
from h5py import File as H5File
from typing import Literal
import shutil

### VIASH START
par = {
    "input": ["resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
              "resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"],
    "output": "foo.h5mu",
    "input_id": ["mouse", "human"],
    "other_axis_mode": "move",
    "output_compression": "gzip"
}
meta = {
    "cpus": 10,
    "resources_dir": "resources_test/"
}
### VIASH END

sys.path.append(meta["resources_dir"])

# START TEMPORARY WORKAROUND compress_h5mu
# reason: resources aren't available when using Nextflow fusion

# from compress_h5mu import compress_h5mu
from h5py import Group, Dataset
from typing import Union
from functools import partial

def compress_h5mu(input_path: Union[str, Path], 
                output_path: Union[str, Path], 
                compression: Union[Literal['gzip'], Literal['lzf']]):
    input_path, output_path = str(input_path), str(output_path)

    def copy_attributes(in_object, out_object):
        for key, value in in_object.attrs.items():
            out_object.attrs[key] = value

    def visit_path(output_h5: H5File,
                   compression: Union[Literal['gzip'], Literal['lzf']], 
                   name: str, object: Union[Group, Dataset]):
            if isinstance(object, Group):
                new_group = output_h5.create_group(name)
                copy_attributes(object, new_group)
            elif isinstance(object, Dataset):
                # Compression only works for non-scalar Dataset objects
                # Scalar objects dont have a shape defined
                if not object.compression and object.shape not in [None, ()]: 
                    new_dataset = output_h5.create_dataset(name, data=object, compression=compression)
                    copy_attributes(object, new_dataset)
                else:
                    output_h5.copy(object, name)
            else:
                raise NotImplementedError(f"Could not copy element {name}, "
                                          f"type has not been implemented yet: {type(object)}")

    with H5File(input_path, 'r') as input_h5, H5File(output_path, 'w', userblock_size=512) as output_h5:
        copy_attributes(input_h5, output_h5)
        input_h5.visititems(partial(visit_path, output_h5, compression))

    with open(input_path, "rb") as input_bytes:
        # Mudata puts metadata like this in the first 512 bytes:
        # MuData (format-version=0.1.0;creator=muon;creator-version=0.2.0)
        # See mudata/_core/io.py, read_h5mu() function
        starting_metadata = input_bytes.read(100)
        # The metadata is padded with extra null bytes up until 512 bytes
        truncate_location = starting_metadata.find(b"\x00")
        starting_metadata = starting_metadata[:truncate_location]
    with open(output_path, "br+") as f:
        nbytes = f.write(starting_metadata)
        f.write(b"\0" * (512 - nbytes))
# END TEMPORARY WORKAROUND compress_h5mu

# START TEMPORARY WORKAROUND setup_logger
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

def indexes_unique(indices: Iterable[pd.Index]) -> bool:
    combined_indices = indices[0].append(indices[1:])
    return combined_indices.is_unique

def check_observations_unique(samples: Iterable[anndata.AnnData]) -> None:
    observation_ids = [sample.obs.index for sample in samples]
    if not indexes_unique(observation_ids):
        raise ValueError("Observations are not unique across samples.")


def nunique(row):
    unique = pd.unique(row)
    unique_without_na = pd.core.dtypes.missing.remove_na_arraylike(unique)
    return len(unique_without_na) > 1

def any_row_contains_duplicate_values(n_processes: int, frame: pd.DataFrame) -> bool:
    """
    Check if any row contains duplicate values, that are not NA.
    """
    numpy_array = frame.to_numpy()
    with Pool(n_processes) as pool:
        is_duplicated = pool.map(nunique, iter(numpy_array))
    return any(is_duplicated)

def concatenate_matrices(n_processes: int, matrices: dict[str, pd.DataFrame]) \
    -> tuple[dict[str, pd.DataFrame], pd.DataFrame | None, dict[str, pd.core.dtypes.dtypes.Dtype]]:
    """
    Merge matrices by combining columns that have the same name.
    Columns that contain conflicting values (e.i. the columns have different values),
    are not merged, but instead moved to a new dataframe.
    """
    column_names = set(column_name for var in matrices.values() for column_name in var)
    logger.debug('Trying to concatenate columns: %s.', ",".join(column_names))
    if not column_names:
        return {}, None
    conflicts, concatenated_matrix = \
        split_conflicts_and_concatenated_columns(n_processes,
                                                 matrices,
                                                 column_names)
    concatenated_matrix = cast_to_writeable_dtype(concatenated_matrix)
    conflicts = {conflict_name: cast_to_writeable_dtype(conflict_df) 
                 for conflict_name, conflict_df in conflicts.items()}
    return conflicts, concatenated_matrix

def get_first_non_na_value_vector(df):
    numpy_arr = df.to_numpy()
    n_rows, n_cols = numpy_arr.shape
    col_index = pd.isna(numpy_arr).argmin(axis=1)
    flat_index = n_cols * np.arange(n_rows) + col_index
    return pd.Series(numpy_arr.ravel()[flat_index], index=df.index, name=df.columns[0])

def split_conflicts_and_concatenated_columns(n_processes: int,
                                             matrices: dict[str, pd.DataFrame],
                                             column_names: Iterable[str]) -> \
                                            tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Retrieve columns with the same name from a list of dataframes which are
    identical across all the frames (ignoring NA values).
    Columns which are not the same are regarded as 'conflicts',
    which are stored in seperate dataframes, one per columns
    with the same name that store conflicting values.
    """
    conflicts = {}
    concatenated_matrix = []
    for column_name in column_names:
        columns = {input_id: var[column_name]
                   for input_id, var in matrices.items()
                   if column_name in var}
        assert columns, "Some columns should have been found."
        concatenated_columns = pd.concat(columns.values(), axis=1, join="outer")
        if any_row_contains_duplicate_values(n_processes, concatenated_columns):
            concatenated_columns.columns = columns.keys() # Use the sample id as column name
            conflicts[f'conflict_{column_name}'] = concatenated_columns
        else:
            unique_values = get_first_non_na_value_vector(concatenated_columns)
            # concatenated_columns.fillna(method='bfill', axis=1).iloc[:, 0]
            concatenated_matrix.append(unique_values)
    if concatenated_matrix:
        concatenated_matrix = pd.concat(concatenated_matrix, join="outer", axis=1)
    else:
        concatenated_matrix = pd.DataFrame()

    return conflicts, concatenated_matrix

def cast_to_writeable_dtype(result: pd.DataFrame) -> pd.DataFrame:
    """
    Cast the dataframe to dtypes that can be written by mudata.
    """    
    # dtype inferral workfs better with np.nan
    result = result.replace({pd.NA: np.nan})

    # MuData supports nullable booleans and ints
    # ie. `IntegerArray` and `BooleanArray`
    result = result.convert_dtypes(infer_objects=True,
                                   convert_integer=True,
                                   convert_string=False,
                                   convert_boolean=True,
                                   convert_floating=False)
    
    # Convert leftover 'object' columns to string
    # However, na values are supported, so convert all values except NA's to string
    object_cols = result.select_dtypes(include='object').columns.values
    for obj_col in object_cols:
        result[obj_col] = result[obj_col].where(result[obj_col].isna(), result[obj_col].astype(str)).astype('category')
    return result

def split_conflicts_modalities(n_processes: int, samples: dict[str, anndata.AnnData], output: anndata.AnnData) \
        -> anndata.AnnData:
    """
    Merge .var and .obs matrices of the anndata objects. Columns are merged
    when the values (excl NA) are the same in each of the matrices.
    Conflicting columns are moved to a separate dataframe (one dataframe for each column,
    containing all the corresponding column from each sample).
    """
    matrices_to_parse = ("var", "obs")
    for matrix_name in matrices_to_parse:
        matrices = {sample_id: getattr(sample, matrix_name) for sample_id, sample in samples.items()}
        conflicts, concatenated_matrix = concatenate_matrices(n_processes, matrices)
        
        # Write the conflicts to the output
        matrix_index = getattr(output, matrix_name).index
        for conflict_name, conflict_data in conflicts.items():
            getattr(output, f"{matrix_name}m")[conflict_name] = conflict_data.reindex(matrix_index)

        # Set other annotation matrices in the output
        setattr(output, matrix_name, pd.DataFrame() if concatenated_matrix is None else concatenated_matrix)

    return output


def concatenate_modality(n_processes: int, mod: str, input_files: Iterable[str | Path], 
                         other_axis_mode: str, input_ids: tuple[str]) -> anndata.AnnData:
    
    concat_modes = {
        "move": None,
    }
    other_axis_mode_to_apply = concat_modes.get(other_axis_mode, other_axis_mode)

    mod_data = {}
    for input_id, input_file in zip(input_ids, input_files):
        try:
            mod_data[input_id] = mu.read_h5ad(input_file, mod=mod)
        except KeyError as e: # Modality does not exist for this sample, skip it
            if f"Unable to open object '{mod}' doesn't exist" not in str(e):
                raise e
            pass
    check_observations_unique(mod_data.values())

    concatenated_data = anndata.concat(mod_data.values(), join='outer', merge=other_axis_mode_to_apply)

    if other_axis_mode == "move":
        concatenated_data = split_conflicts_modalities(n_processes, mod_data, concatenated_data)
    
    return concatenated_data

def concatenate_modalities(n_processes: int, modalities: list[str], input_files: Path | str,
                           other_axis_mode: str, output_file: Path | str,
                           compression: Literal['gzip'] | Literal['lzf'],
                           input_ids: tuple[str] | None = None) -> None:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info('Concatenating samples.')
    output_file, input_files = Path(output_file), [Path(input_file) for input_file in input_files]
    output_file_uncompressed = output_file.with_name(output_file.stem + "_uncompressed.h5mu")
    output_file_uncompressed.touch()
    # Create empty mudata file
    mdata = mu.MuData({modality: anndata.AnnData() for modality in modalities})
    mdata.write(output_file_uncompressed, compression=compression)

    for mod_name in modalities:
        new_mod = concatenate_modality(n_processes, mod_name, 
                                       input_files, other_axis_mode, 
                                       input_ids)
        logger.info("Writing out modality '%s' to '%s' with compression '%s'.",
                    mod_name, output_file_uncompressed, compression)
        mu.write_h5ad(output_file_uncompressed, data=new_mod, mod=mod_name)
    
    if compression:
        compress_h5mu(output_file_uncompressed, output_file, compression=compression)
        output_file_uncompressed.unlink()
    else:
        shutil.move(output_file_uncompressed, output_file)

    logger.info("Concatenation successful.")

def main() -> None:
    # Get a list of all possible modalities
    mods = set()
    for path in par["input"]:
        try:
            with H5File(path, 'r') as f_root:
                mods = mods | set(f_root["mod"].keys())
        except OSError:
            raise OSError(f"Failed to load {path}. Is it a valid h5 file?")

    input_ids = None
    if par["input_id"]:
        input_ids: tuple[str] = tuple(i.strip() for i in par["input_id"])
        if len(input_ids) != len(par["input"]):
            raise ValueError("The number of sample names must match the number of sample files.")

        if len(set(input_ids)) != len(input_ids):
            raise ValueError("The sample names should be unique.")

    logger.info("\nConcatenating data from paths:\n\t%s",
                "\n\t".join(par["input"]))

    if par["other_axis_mode"] == "move" and not input_ids:
        raise ValueError("--mode 'move' requires --input_ids.")

    n_processes = meta["cpus"] if meta["cpus"] else 1
    concatenate_modalities(n_processes,
                           list(mods),
                           par["input"],
                           par["other_axis_mode"],
                           par["output"],
                           par["output_compression"],
                           input_ids=input_ids)


if __name__ == "__main__":
    main()
