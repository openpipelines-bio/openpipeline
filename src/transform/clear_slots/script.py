import sys

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
from setup_logger import setup_logger


logger = setup_logger()


def main(par):
    logger.info("clear_slots scaffold placeholder")


if __name__ == "__main__":
    sys.exit(main(par))
