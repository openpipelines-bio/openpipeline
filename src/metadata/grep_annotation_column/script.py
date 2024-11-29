import sys
import mudata as mu
from pathlib import Path
from operator import attrgetter
from pandas import Series
import scipy as sc
import re
import numpy as np


### VIASH START
par = {
    "input": "./resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    "modality": "rna",
    "matrix": "var",
    "input_column": "gene_symbol",
    "regex_pattern": "^[mM][tT]-",
    "output": "foo.h5mu",
    "input_id": "mouse",
    "output_match_column": "test",
    "output_fraction_column": "fraction_test",
    "output_compression": "gzip"
}
### VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

def describe_array(arr, msg):
    # Note: sc.stats returns a DescribeResult NamedTuple. For NamedTuples,
    # the _asdict method is public facing even though it starts with an underscore.
    description = sc.stats.describe(arr)._asdict()
    logger.info("%s:\nshape: %s\nmean: %s\nnobs: %s\n"
                "variance: %s\nmin: %s\nmax: %s\ncontains na: %s\ndtype: %s\ncontains 0: %s",
                msg, arr.shape, description["mean"], description["nobs"],
                description["variance"], description["minmax"][0],
                description["minmax"][1], np.isnan(arr).any(), arr.dtype,
                (arr == 0).any())


def main(par):
    input_file, output_file, mod_name = Path(par["input"]), Path(par["output"]), par['modality']
    logger.info(f"Compiling regular expression '{par['regex_pattern']}'.")
    try:
        compiled_regex = re.compile(par["regex_pattern"])
    except (TypeError, re.error) as e:
        raise ValueError(f"{par['regex_pattern']} is not a valid regular expression pattern.") from e
    else:
        if compiled_regex.groups:
            raise NotImplementedError("Using match groups is not supported by this component.")
    logger.info('Reading input file %s, modality %s.', input_file, mod_name)

    mudata = mu.read_h5mu(input_file)
    modality_data = mudata[mod_name]
    logger.info("Reading input file done.")
    logger.info("Using annotation dataframe '%s'.", par["matrix"])
    annotation_matrix = getattr(modality_data, par['matrix'])
    default_column = {
        "var": attrgetter("var_names"),
        "obs": attrgetter("obs_names")
    }
    if par["input_column"]:
        logger.info("Input column '%s' was specified.", par["input_column"])
        try:
            annotation_column = annotation_matrix[par["input_column"]]
        except KeyError as e:
            raise ValueError(f"Column {par['input_column']} could not be found for modality "
                            f"{par['modality']}. Available columns:"
                            f" {','.join(annotation_matrix.columns.to_list())}") from e
    else:
        logger.info(f"No input column specified, using '.{par['matrix']}_names'")
        annotation_column = default_column[par['matrix']](modality_data).to_series()
    logger.info("Applying regex search.")
    grep_result = annotation_column.str.contains(par["regex_pattern"], regex=True)
    logger.info("Search results: %s", grep_result.value_counts())

    other_axis_attribute = {
        "var": "obs",
        "obs": "var"
    }
    if par['output_fraction_column']:
        logger.info("Enabled writing the fraction of values that matches to the pattern.")
        input_layer = modality_data.X if not par["input_layer"] else modality_data.layers[par["input_layer"]]
        totals = np.ravel(input_layer.sum(axis=1))
        describe_array(totals, "Summary of total counts for layer")
        counts_for_matches = np.ravel(input_layer[:, grep_result].sum(axis=1))
        describe_array(counts_for_matches, "Summary of counts matching grep")
        with np.errstate(all='raise'):
            pct_matching = np.divide(counts_for_matches, totals,
                                     out=np.zeros_like(totals, dtype=np.float64),
                                     where=(~np.isclose(totals, np.zeros_like(totals))))
        logger.info("Testing wether or not fractions data contains NA.")
        assert ~np.isnan(pct_matching).any(), "Fractions should not contain NA."
        logger.info("Fraction statistics: \n%s", Series(pct_matching).describe())
        pct_matching = np.where(np.isclose(pct_matching, 0, atol=1e-6), 0, pct_matching)
        pct_matching = np.where(np.isclose(pct_matching, 1, atol=1e-6), 1, pct_matching)
        assert (np.logical_and(pct_matching >= 0, pct_matching <= 1)).all(), \
                "Fractions are not within bounds, please report this as a bug"
        output_matrix = other_axis_attribute[par['matrix']]
        logger.info("Writing fractions to matrix '%s', column '%s'",
                    output_matrix, par['output_fraction_column'])
        getattr(modality_data, output_matrix)[par['output_fraction_column']] = pct_matching
    logger.info("Adding values that matched the pattern to '%s', column '%s'",
                par["matrix"], par["output_match_column"])
    getattr(modality_data, par['matrix'])[par["output_match_column"]] = grep_result
    logger.info("Writing out data to '%s' with compression '%s'.",
                output_file, par["output_compression"])
    mudata.write(output_file, compression=par["output_compression"])

if __name__ == "__main__":
    main(par)