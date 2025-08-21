import sys
import mudata
import tiledbsoma
import tiledbsoma.io
import pandas as pd
import os
import json

## VIASH START
par = {}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

tiledbsoma.logging.info()
print(os.environ, file=sys.stderr, flush=True)


def main(par):
    logger.info(f"Component {meta['name']} started.")
    par["input_uri"] = par["input_uri"].rstrip("/")
    if not par["obsm_input"]:
        raise ValueError("Please provide at least one .obsm key.")
    logger.info(
        "Opening mudata file '%s', modality '%s'.", par["input_mudata"], par["modality"]
    )
    modality_data = mudata.read_h5ad(par["input_mudata"], mod=par["modality"])
    logger.info(
        "Done reading modality. Looking at .obsm for keys: '%s'",
        ",".join(par["obsm_input"]),
    )
    try:
        keys_to_transfer = {
            obsm_key: modality_data.obsm[obsm_key] for obsm_key in par["obsm_input"]
        }
    except KeyError as e:
        raise KeyError("Not all .obsm keys were found in the input!") from e

    logger.info("Done getting .obsm keys.")
    tiledb_config = {
        "vfs.s3.no_sign_request": "false",
        "vfs.s3.region": par["s3_region"],
    }
    optional_config = {
        "vfs.s3.endpoint_override": par["endpoint"],
    }
    for config_setting, config_val in optional_config.items():
        if config_val is not None:
            tiledb_config[config_setting] = config_val
    logger.info("Using the following config to connect: %s", tiledb_config)
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)

    logger.info(
        "Trying to access '%s' in region '%s'", par["input_uri"], par["s3_region"]
    )
    with tiledbsoma.open(
        par["input_uri"], mode="w", context=context
    ) as open_experiment:
        logger.info("Connection established.")
        logger.info("Looking for measurement %s", par["output_modality"])
        measurement = open_experiment.ms[par["output_modality"]]
        logger.info("Checking if keys do not already exist.")
        existing_keys = measurement.obsm.keys()
        overlap = set(existing_keys).intersection(set(keys_to_transfer))
        if overlap:
            raise ValueError(
                f"The following keys already exist in the database: {','.join(overlap)}."
            )
        logger.info("Adding keys to database.")
        for key, obsm_val in keys_to_transfer.items():
            logger.info("Adding .obsm key '%s', of class '%s'", key, type(obsm_val))
            index_as_json = None
            if isinstance(obsm_val, pd.DataFrame):
                # tileDB does not allow column indices to be saved directly
                # So need to add those as JSON metadata
                index_to_write = obsm_val.columns.to_list()
                if not isinstance(index_to_write, pd.RangeIndex):
                    index_as_json = json.dumps(index_to_write)
                obsm_val = obsm_val.to_numpy()

            tiledbsoma.io.add_matrix_to_collection(
                open_experiment,
                measurement_name=par["output_modality"],
                ingest_mode="write",
                collection_name="obsm",
                matrix_name=key,
                matrix_data=obsm_val,
                context=context,
            )
            if index_as_json:
                uri = f"{par['input_uri']}/ms/{par['output_modality']}/obsm/{key}"
                with tiledbsoma.open(
                    uri=uri, mode="w", context=context
                ) as open_obsm_array:
                    open_obsm_array.metadata["column_index"] = index_as_json

    logger.info("Finished!")


if __name__ == "__main__":
    main(par)
