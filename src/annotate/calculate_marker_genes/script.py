import sys

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "input_layer": None,
    "sel_clustering": "leiden_0.5",
    "output": "output.h5mu",
    "output_markers": "markers.csv",
    "output_compression": "gzip",
    "method": "t-test",
    "corr_method": "benjamini-hochberg",
    "tie_correct": False,
    "key_added": "rank_genes_groups",
    "filter_results": False,
    "min_in_group_fraction": 0.25,
    "max_out_group_fraction": 0.5,
    "min_fold_change": 1,
    "compare_abs": False,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger


logger = setup_logger()


def main(par):
    logger.info("calculate_marker_genes scaffold placeholder")


if __name__ == "__main__":
    sys.exit(main(par))