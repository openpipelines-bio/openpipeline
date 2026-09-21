import mudata as mu
from compress_h5mu import write_h5ad_to_h5mu_with_compression


def read_modality(input_file, modality, logger):
    logger.info("Reading modality %s from %s", modality, input_file)
    dat = mu.read_h5ad(input_file, mod=modality)
    assert dat.var_names.is_unique, (
        "The var_names of the input modality must be unique."
    )
    return dat


def write_modality(dat, output_file, input_file, modality, output_compression, logger):
    logger.info(
        "Writing to file to %s with compression %s",
        output_file,
        output_compression,
    )
    write_h5ad_to_h5mu_with_compression(
        output_file, input_file, modality, dat, output_compression
    )
