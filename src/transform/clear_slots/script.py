import sys
from pathlib import Path

from mudata import read_h5ad

## VIASH START
par = {
    "input": "input.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "slots": ["obsm", "varm"],
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import write_h5ad_to_h5mu_with_compression
from setup_logger import setup_logger


logger = setup_logger()


def main(par):
    input_file = Path(par["input"])
    output_file = Path(par["output"])
    mod_name = par["modality"]

    logger.info("Reading input file %s, modality %s.", input_file, mod_name)
    mod = read_h5ad(input_file, mod=mod_name)

    for slot in par["slots"]:
        if not hasattr(mod, slot):
            raise ValueError(f"Slot '{slot}' is not supported for modality {mod_name}.")

        logger.info("Clearing %s from modality %s.", slot, mod_name)
        getattr(mod, slot).clear()

    logger.info("Writing output to %s.", output_file)
    write_h5ad_to_h5mu_with_compression(
        output_file, input_file, mod_name, mod, par["output_compression"]
    )

    logger.info("Done.")


if __name__ == "__main__":
    sys.exit(main(par))
