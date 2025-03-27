from __future__ import annotations
import sys
import mudata as md
import pandas as pd
import numpy as np


### VIASH START
par = {
    "input": [
        "./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu",
        "./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu",
    ],
    "output": "foo.h5mu",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main():
    logger.info("Reading input files %s", ",".join(par["input"]))
    input_samples = [md.read_h5mu(path) for path in par["input"]]

    logger.info("Merging into single object.")
    sample_modalities = {}
    for input_sample in input_samples:
        for mod_name, mod_data in input_sample.mod.items():
            if mod_name in sample_modalities:
                raise ValueError(
                    f"Modality '{mod_name}' was found in more than 1 sample."
                )
            sample_modalities[mod_name] = mod_data

    merged = md.MuData(sample_modalities)
    merged.update()
    for df_attr in ("var", "obs"):
        df = getattr(merged, df_attr)
        df = df.replace({pd.NA: np.nan}, inplace=False)

        # MuData supports nullable booleans and ints
        # ie. `IntegerArray` and `BooleanArray`
        df = df.convert_dtypes(
            infer_objects=True,
            convert_integer=True,
            convert_string=False,
            convert_boolean=True,
            convert_floating=False,
        )

        # pd.convert_dtypes does not convert already existing nullable dtypes
        # to their native numpy dtypes.
        # At this point, integer and boolean columns use a nullable dtype; and string columns
        # might use 'object' or pd.StringDtype. This is OK.
        # However, for example floating number columns might still use a nullable dtype
        # from before the call to convert_dtypes. So we need to explicitly cast.
        col_dtypes = df.select_dtypes(
            exclude=["integer", "string", "object", "boolean"]
        ).dtypes.to_dict()
        numpy_dtypes = {
            col_name: np.dtype(col_dtype.kind)
            for col_name, col_dtype in col_dtypes.items()
            if pd.api.types.is_extension_array_dtype(col_dtype)
        }
        df = df.astype(numpy_dtypes)

        # Convert leftover 'object' columns to string
        object_cols = df.select_dtypes(include="object").columns.values
        for obj_col in object_cols:
            df[obj_col] = df[obj_col].astype(str).astype("category")

        setattr(merged, df_attr, df)

    merged.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished")


if __name__ == "__main__":
    main()
