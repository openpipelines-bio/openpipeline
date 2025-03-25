from __future__ import annotations
import sys
import mudata as md
import pandas as pd
import numpy as np


### VIASH START
par = {
    "input": [
        "merged.align_query_reference.output_query_prot.h5mu",
        "merged.knn.output.h5mu",
    ],
    "output": "foo.h5mu",
    "output_compression": None
}
meta = {
    "resources_dir": "src/utils"
}
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
            convert_floating=True,
        )

        # Convert leftover 'object' columns to string
        object_cols = df.select_dtypes(include="object").columns.values
        for obj_col in object_cols:
            df[obj_col] = df[obj_col].astype(str).astype("category")
            
        # Nullable float columns are not supported by anndata/mudata
        float64_cols = df.select_dtypes(include=pd.Float64Dtype())
        for float64_col in float64_cols:
            df[float64_col] = df[float64_col].astype(np.float64)
        
        float32_cols = df.select_dtypes(include=pd.Float32Dtype())
        for float32_col in float32_cols:
            df[float32_col] = df[float32_col].astype(np.float32)
            
        setattr(merged, df_attr, df)

    merged.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished")


if __name__ == "__main__":
    main()
