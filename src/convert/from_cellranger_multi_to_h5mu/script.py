from pathlib import Path
import logging
from sys import stdout
import scanpy
import pandas as pd
import mudata
from scirpy.io import read_10x_vdj
from collections import defaultdict
## VIASH START
par = {
    "input": "resources_test/10x_5k_anticmv/processed/10x_5k_anticmv.cellranger_multi.output.output",
    "output": "foo.h5mu",
    "uns_metrics": "metrics_cellranger",
    "output_compression": "gzip"
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

POSSIBLE_LIBRARY_TYPES = ('vdj_t', 'vdj_b', 'vdj_t_gd', 'count')


def gather_input_data(dir: Path):
    if not dir.is_dir():
        raise ValueError("Specified input is not a directory.")
    folder_contents = list(dir.iterdir())
    config = dir / 'config.csv'
    if config not in folder_contents:
        logger.warning('Config.csv not found in input directory, this folder might not be a valid cellranger multi output.')

    required_subfolders = [dir / subfolder_name for subfolder_name in ('multi', 'per_sample_outs')]
    found_input = {key_: None for key_ in POSSIBLE_LIBRARY_TYPES + ('metrics_summary',)}
    for required_subfolder in required_subfolders:
        if not required_subfolder in folder_contents:
            raise ValueError(f"Input folder must contain the subfolder {required_subfolder} please make "
                               "sure that the specified input folder is a valid cellranger multi output.")

    multi_dir = dir / 'multi'
    for library_type in multi_dir.iterdir():
        if not library_type.is_dir():
            logger.warning("%s is not a directory. Contents of the multi folder "
                           "must be directories to be recognized as valid input data",
                           library_type)
            continue
        if library_type.name not in POSSIBLE_LIBRARY_TYPES:
            raise ValueError(f"Contents of the 'multi' folder must be found one of the following: {','.join(POSSIBLE_LIBRARY_TYPES)}.")

        found_input[library_type.name] = library_type

    per_sample_outs_dir = dir / 'per_sample_outs'
    metric_summaries = list(per_sample_outs_dir.glob('*/metrics_summary.csv'))
    if len(metric_summaries) > 1:
        raise ValueError("Found more than one 'metrics_summery.csv' file. "
                         "This component currently only supports parsing cellranger multi output for one sample.")
    found_input["metrics_summary"] = metric_summaries[0]
    return found_input


def process_counts(counts_folder: Path):
    counts_matrix_file = counts_folder / "raw_feature_bc_matrix.h5"
    logger.info("Reading %s.", counts_matrix_file)
    adata = scanpy.read_10x_h5(counts_matrix_file, gex_only=False)

    # set the gene ids as var_names
    logger.info("Renaming var columns")
    adata.var = adata.var\
        .rename_axis("gene_symbol")\
        .reset_index()\
        .set_index("gene_ids")

    # generate output
    logger.info("Convert to mudata")

    def modality_name_factory(library_type):
        return ("".join(library_type.replace("-", "_").split())).lower()

    feature_types = defaultdict(modality_name_factory, {
            "Gene Expression": "rna",
            "Peaks": "atac",
            "Antibody Capture": "prot",
            "VDJ": "vdj",
            "VDJ-T": "vdj_t",
            "VDJ-B": "vdj_b",
            "CRISPR Guide Capture": "gdo",
            "Multiplexing Capture": "hto"
        })
    return mudata.MuData(adata, feature_types=feature_types)

def process_metrics_summary(mudata: mudata.MuData, metrics_file: Path):
    def read_percentage(val):
        try:
            return float(val.strip('%')) / 100
        except (AttributeError, ValueError):
            return val

    metrics_summary = pd.read_csv(metrics_file,
                                  decimal=".",
                                  quotechar='"',
                                  thousands=",").applymap(read_percentage)

    mudata.uns[par["uns_metrics"]] = metrics_summary
    for colname, coldata in metrics_summary.items():
        try:
            new_column = coldata.astype(str, copy=True).astype({colname: "category"})
            metrics_summary[colname] = new_column
        except (ValueError, TypeError):
            logger.warning(f"Could not store column {colname} from metrics.")
            pass
    return mudata

def process_vdj(mudata: mudata.MuData, vdj_folder_path: Path):
    # https://scverse.org/scirpy/latest/generated/scirpy.io.read_10x_vdj.html#scirpy-io-read-10x-vdj
    # According to docs, using the json is preferred as this file includes intron info.
    all_config_json_file = vdj_folder_path / "all_contig_annotations.json"
    vdj_anndata = read_10x_vdj(all_config_json_file)
    vdj_type = vdj_folder_path.name
    mudata.mod[vdj_type] = vdj_anndata
    return mudata

def get_modalities(input_data):
    dispatcher = {
        'vdj_t': process_vdj,
        'vdj_b': process_vdj,
        'vdj_t_gd': process_vdj,
        'metrics_summary': process_metrics_summary
    }
    mudata_file = process_counts(input_data['count'])
    for modality_name, modality_data_path in input_data.items():
        if modality_name == "count" or not modality_data_path:
            continue
        try:
            parser_function = dispatcher[modality_name]
        except KeyError as e:
            raise ValueError("This component does not support the "
                             f"parsing of the modality '{modality_name}' yet.") from e
        mudata_file = parser_function(mudata_file, modality_data_path)
    return mudata_file

def main():
    cellranger_multi_dir = Path(par["input"])
    input_data = gather_input_data(cellranger_multi_dir)
    result = get_modalities(input_data)
    logger.info("Writing %s", par["output"])
    result.write_h5mu(par["output"], compression=par["output_compression"])

if __name__ == "__main__":
    main()