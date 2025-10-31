import sys
import mudata
import tiledbsoma
import tiledbsoma.io
import pandas as pd
from pathlib import Path

## VIASH START
par = {
    "input_uri": "s3://openpipelines-data/tiledb/pbmc_1k_protein_v3_mms/",
    "input_mudata": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "obs_input": ["test_slot"],
    "modality": "rna",
    "input_dir": None,
    "s3_region": "eu-west-3",
    "endpoint": None,
    "output_tiledb": "./output",
    "output_modality": "rna",
    "s3_no_sign_request": True,
    "obs_index_name_input": "cell_id",
}
meta = {"resources_dir": "src/utils", "name": "move_mudata_obs_to_tiledb"}

test_path = "./mudata_for_testing.h5mu"
test_mudata = mudata.read_h5mu(par["input_mudata"])
test_mudata["rna"].obs["test_slot"] = (
    test_mudata["rna"].obs["filter_with_counts"].copy()
)
test_mudata.write(test_path)
par["input_mudata"] = test_path
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

tiledbsoma.logging.info()


def download_s3_dir(input_uri, output_dir_path):
    logger.info("Requested to output to a directory. Downloading database...")
    output_dir_path.mkdir(parents=True, exist_ok=True)
    import boto3
    from awscli.customizations.s3.utils import split_s3_bucket_key

    bucket, key = split_s3_bucket_key(input_uri)
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


def main(par):
    logger.info(f"Component {meta['name']} started.")
    if par["input_uri"]:
        par["input_uri"] = par["input_uri"].rstrip("/")
    if par["input_uri"] and par["input_dir"]:
        raise ValueError("Cannot provide both 'input_uri' and 'input_dir'.")
    if not par["input_uri"] and not par["input_dir"]:
        raise ValueError("Must provide either 'input_uri' or 'input_dir'")
    if not par["obs_input"]:
        raise ValueError("Please provide at least one .obs column.")
    logger.info(
        "Opening mudata file '%s', modality '%s'.", par["input_mudata"], par["modality"]
    )
    modality_data = mudata.read_h5ad(par["input_mudata"], mod=par["modality"])
    logger.info(
        "Done reading modality. Looking at .obs for keys: '%s'",
        ",".join(par["obs_input"]),
    )
    try:
        for obs_key in par["obs_input"]:
            modality_data.obs[obs_key]
    except KeyError as e:
        raise KeyError("Not all .obs keys were found in the input!") from e

    logger.info("Done getting .obs keys.")
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

    if par["output_tiledb"]:
        output_dir_path = Path(par["output_tiledb"])
        if par["input_dir"]:
            import shutil

            shutil.copytree(par["input_dir"], output_dir_path, dirs_exist_ok=True)
        else:
            download_s3_dir(par["input_uri"], output_dir_path)
        logger.info("Setting input to '%s'", output_dir_path)
        par["input_uri"] = f"file://{output_dir_path.resolve()}"
        logger.info("Overwriting TileDB config because S3 connection is not required.")
        tiledb_config = {}

    logger.info("Trying to access '%s'", par["input_uri"])
    logger.info("Fetching .obs")
    with tiledbsoma.open(
        par["input_uri"],
        mode="r",
        context=tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config),
    ) as open_experiment:
        logger.info("Connection established.")
        obs_df = open_experiment.obs.read().concat().to_pandas()
        logger.info("Done downloading .obs from databse.")
    logger.info("Adding obs columns to fetched .obs dataframe.")
    overlapping_obs = set(par["obs_input"]).intersection(set(obs_df.columns.to_list()))
    if overlapping_obs:
        raise ValueError(
            f"The following keys already exist in the database: {','.join(overlapping_obs)}."
        )

    columns_to_add = modality_data.obs[par["obs_input"]]
    new_obs = pd.merge(
        obs_df,
        columns_to_add,
        left_on=par["obs_index_name_input"],
        right_index=True,
        how="right",
    )
    logger.info(
        "Writing obs back to database. Connection to %s with config %s",
        par["input_uri"],
        tiledb_config,
    )
    with tiledbsoma.open(
        par["input_uri"],
        mode="w",
        context=tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config),
    ) as open_experiment:
        tiledbsoma.io.update_obs(
            open_experiment,
            new_data=new_obs,
            default_index_name=par["obs_index_name_input"],
            context=tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config),
        )
    logger.info("Finished!")


if __name__ == "__main__":
    main(par)
