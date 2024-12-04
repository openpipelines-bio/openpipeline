import sys
from mudata import read_h5ad, write_h5ad
import shutil
from pathlib import Path

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": ["log_normalized"],
    "missing_ok": False,
    "output_compression": "lzf",
}
meta = {"name": "delete_layer", "resources_dir": "resources_test"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import compress_h5mu

logger = setup_logger()


def main():
    input_file, output_file, mod_name = (
        Path(par["input"]),
        Path(par["output"]),
        par["modality"],
    )

    logger.info("Reading input file %s, modality %s.", input_file, mod_name)
    mod = read_h5ad(input_file, mod=mod_name)
    for layer in par["layer"]:
        if layer not in mod.layers:
            if par["missing_ok"]:
                continue
            raise ValueError(f"Layer '{layer}' is not present in modality {mod_name}.")
        logger.info("Deleting layer %s from modality %s.", layer, mod_name)
        del mod.layers[layer]

    logger.info("Writing output to %s.", par["output"])
    output_file_uncompressed = (
        output_file.with_name(output_file.stem + "_uncompressed.h5mu")
        if par["output_compression"]
        else output_file
    )
    shutil.copyfile(par["input"], output_file_uncompressed)
    write_h5ad(filename=output_file_uncompressed, mod=mod_name, data=mod)
    if par["output_compression"]:
        compress_h5mu(
            output_file_uncompressed, output_file, compression=par["output_compression"]
        )
        output_file_uncompressed.unlink()

    logger.info("Finished.")


if __name__ == "__main__":
    main()
