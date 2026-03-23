from __future__ import annotations
from functools import partial
import sys
import mudata as md
import anndata as ad
from pathlib import Path
import pandas as pd


### VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo/",
    "output_types": "foo_types.csv",
    "output_compression": "gzip",
}
meta = {"resources_dir": "./src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from mudata_opener import mudata_opener

logger = setup_logger()


def main() -> None:
    output_dir = Path(par["output"])
    logger.info("Creating output directory '%s' if it does not exist", output_dir)
    if not output_dir.is_dir():
        logger.info("Creating %s", output_dir)
        output_dir.mkdir(parents=True)

    input_file = Path(par["input"])
    logger.info("Checking which modalities exist for '%s'", par["input"])
    with mudata_opener(input_file, mode="r") as (open_mudata, input_is_zarr):
        logger.info(
            "Openened %s in %s format.", par["input"], "zarr" if input_is_zarr else "h5"
        )

        modalities = list(open_mudata["mod"].keys())
        logger.info("Found the following modalities:\n%s", "\n".join(modalities))

        logger.info("Creating output types CSV.")
        output_extension = "zarr" if input_is_zarr else "h5mu"
        names = {
            mod_name: f"{input_file.stem}_{mod_name}.{output_extension}"
            for mod_name in modalities
        }
        output_files = list(names.values())
        logger.info(
            "Will be creating the following output .%s files:\n%s",
            output_extension,
            "\n".join(output_files),
        )
        df = pd.DataFrame({"name": modalities, "filename": output_files})
        logger.info("Writing output_types CSV file to '%s'.", par["output_types"])
        df.to_csv(par["output_types"], index=False)

        logger.info("Splitting input file into unimodal output files.")
        for mod_name in modalities:
            logger.info("Processing modality '%s'", mod_name)
            elem_key = f"/mod/{mod_name}"
            elem = open_mudata[elem_key]
            logger.info("Reading %s", elem_key)
            new_ad = ad.io.read_elem(elem)
            logger.info("Creating MuData object.")
            new_sample = md.MuData({mod_name: new_ad})
            logger.info(
                "Writing to '%s', with compression '%s'",
                names[mod_name],
                par["output_compression"],
            )
            writer = (
                partial(md.MuData.write_zarr, zarr_format=3)
                if input_is_zarr
                else partial(
                    md.MuData.write_h5mu, compression=par["output_compression"]
                )
            )
            writer(new_sample, output_dir / names[mod_name])
        logger.info("Done writing output file.")
    logger.info("Finished")


if __name__ == "__main__":
    main()
