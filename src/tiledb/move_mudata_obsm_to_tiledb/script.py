import sys
import mudata
import tiledbsoma
import tiledbsoma.io
import pandas as pd
import json
from pathlib import Path

## VIASH START
par = {
    "input_uri": "s3://openpipelines-data/tiledb/pbmc_1k_protein_v3_mms/",
    "input_mudata": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "obsm_input": ["test_slot"],
    "modality": "rna",
    "s3_region": "eu-west-3",
    "endpoint": None,
    "output_tiledb": "./output",
    "output_modality": "rna",
}
meta = {"resources_dir": "src/utils", "name": "move_mudata_obsm_to_tiledb"}

test_path = "./mudata_for_testing.h5mu"
test_mudata = mudata.read_h5mu(par["input_mudata"])
test_mudata["rna"].obsm["test_slot"] = test_mudata["rna"].obsm["X_pca"].copy()
test_mudata.write(test_path)
par["input_mudata"] = test_path
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

tiledbsoma.logging.info()


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
    optional_config = {
        "vfs.s3.region": par["s3_region"],
        "vfs.s3.endpoint_override": par["endpoint"],
        "vfs.s3.no_sign_request": par["s3_no_sign_request"],
    }
    tiledb_config = {}
    for config_setting, config_val in optional_config.items():
        if config_val is not None:
            tiledb_config[config_setting] = config_val
    logger.info("Using the following config to connect to S3: %s", tiledb_config)
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)

    if par["output_tiledb"]:
        logger.info("Requested to output to a directory. Downloading database...")
        output_dir_path = Path(par["output_tiledb"])
        output_dir_path.mkdir(parents=True, exist_ok=True)
        import boto3
        from awscli.customizations.s3.utils import split_s3_bucket_key

        bucket, key = split_s3_bucket_key(par["input_uri"])
        connection_args = {
            "endpoint_url": par["endpoint"],
            "region_name": par["s3_region"],
        }
        if par["s3_no_sign_request"]:
            import botocore

            connection_args["config"] = botocore.config.Config(
                signature_version=botocore.UNSIGNED
            )

        client = boto3.resource("s3", **connection_args)
        bucket = client.Bucket(bucket)
        for i, s3_obj in enumerate(bucket.objects.filter(Prefix=key)):
            output_path = output_dir_path / s3_obj.key.removeprefix(key).lstrip("/")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            bucket.download_file(s3_obj.key, output_path)
            print(f"Downloaded {i} files.", file=sys.stdout, flush=True, end="\r")
        logger.info("Download completed!")
        logger.info("Setting input to %s", par["input_uri"])
        par["input_uri"] = f"file://{output_dir_path.resolve()}"
        logger.info("Overwriting TileDB config because S3 connection is not required.")
        context = tiledbsoma.SOMATileDBContext()

    logger.info("Trying to access '%s'", par["input_uri"])
    with tiledbsoma.open(
        par["input_uri"], mode="w", context=context
    ) as open_experiment:
        logger.info("Connection established.")
        logger.info("Looking for measurement %s", par["output_modality"])
        measurement = open_experiment.ms[par["output_modality"]]
        logger.info("Checking if keys do not already exist.")
        existing_keys = measurement.obsm.keys()
        logger.info("Existing keys: %s", ",".join(existing_keys))
        overlap = set(existing_keys).intersection(set(keys_to_transfer.keys()))
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
            # Allow for more than one obsm slot to be transferred
            open_experiment = open_experiment.reopen(mode="w")
            if index_as_json:
                uri = f"{par['input_uri']}/ms/{par['output_modality']}/obsm/{key}"
                with tiledbsoma.open(
                    uri=uri, mode="w", context=context
                ) as open_obsm_array:
                    open_obsm_array.metadata["column_index"] = index_as_json

    logger.info("Finished!")


if __name__ == "__main__":
    main(par)
