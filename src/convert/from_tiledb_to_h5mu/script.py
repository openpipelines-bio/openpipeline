import tiledbsoma
import tiledbsoma.io
import sys
import mudata
import json
import pandas as pd

## VIASH START
par = {
    "input_uri": "s3://openpipelines-data/tiledb/pbmc_1k_protein_v3_mms/",
    "input_modality": "rna",
    "output": "output.h5mu",
    "output_modality": "rna",
    "s3_region": "eu-west-3",
    "endpoint": None,
    "input_layers": ["raw", "log_normalized"],
    "output_compression": "gzip",
}

meta = {
    "resources_dir": "src/utils",
    "name": "from_tiledb_to_h5mu",
}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
tiledbsoma.logging.info()


def _log_arguments(function_obj, arg_dict):
    """
    Format a dictionairy of arguments into a string that is put into the script logs.
    """
    args_str = [f"\t{param}: {param_val}\n" for param, param_val in arg_dict.items()]
    logger.info(
        "Calling %s with arguments:\n%s",
        function_obj.__name__,
        "".join(args_str).rstrip(),
    )


def main(par):
    logger.info("Component %s started", meta["name"])
    tiledb_config = {
        "vfs.s3.no_sign_request": "false",
    }
    optional_config = {
        "vfs.s3.region": par["s3_region"],
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
        par["input_uri"], mode="r", context=context
    ) as open_experiment:
        logger.info("Connection successful")
        to_anndata_args = {
            "experiment": open_experiment,
            "measurement_name": par["input_modality"],
            "extra_X_layer_names": par["input_layers"],
        }
        func_to_call = tiledbsoma.io.to_anndata
        _log_arguments(func_to_call, to_anndata_args)
        output_modality = func_to_call(**to_anndata_args)
        logger.info("Output anndata was:\n%s", output_modality)
        logger.info("Adding column indices to varm and obsm items.")
        for multimodal_attribute in ("varm", "obsm"):
            multimodal_dict = getattr(output_modality, multimodal_attribute)
            for key_ in list(multimodal_dict.keys()):
                metadata = getattr(
                    open_experiment.ms[par["input_modality"]], multimodal_attribute
                )[key_].metadata
                if "column_index" in metadata:
                    logger.info(
                        "Found index for item %s in %s", key_, multimodal_attribute
                    )
                    index_list = json.loads(metadata["column_index"])
                    assert isinstance(index_list, list)
                    multimodal_dict[key_] = pd.DataFrame(
                        multimodal_dict[key_],
                        columns=index_list,
                        index=multimodal_dict.dim_names,
                    )
                logger.info(
                    "No column index found for item %s in %s",
                    key_,
                    multimodal_attribute,
                )

    logger.info("Converting to MuData")
    output_mudata = mudata.MuData({par["output_modality"]: output_modality})
    logger.info(
        "Writing output MuData to %s with compression %s",
        par["output"],
        par["output_compression"],
    )
    output_mudata.write(par["output"], compression=par["output_compression"])
    logger.info("Finished!")


if __name__ == "__main__":
    main(par)
