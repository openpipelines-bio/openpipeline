import sys
import tiledbsoma

## VIASH START
par = {
    "input_uri": "s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/",
    "s3_region": "us-west-2",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main(par):
    logger.info("Component started")

    tiledb_config = {
        "vfs.s3.no_sign_request": "false",
        "vfs.s3.region": par["s3_region"],
    }
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)
    logger.info(
        "Trying to access '%s' in region '%s'", par["input_uri"], par["s3_region"]
    )
    with tiledbsoma.open(par["input_uri"], "w", context=context) as _:
        logger.info("Connection successful!")
    logger.info("Component finished")


if __name__ == "__main__":
    main(par)
