from pathlib import Path
import sys
import scanpy
import pandas as pd
import mudata
import numpy as np
from scirpy.io import read_10x_vdj
from collections import defaultdict
from functools import partial


## VIASH START
par = {
    "input": "resources_test/10x_5k_beam/processed/10x_5k_beam.cellranger_multi.output",
    "output": "foo.h5mu",
    "uns_metrics": "metrics_cellranger",
    "output_compression": "gzip"
}
meta = {
    "resources_dir": "."
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

POSSIBLE_LIBRARY_TYPES = ('vdj_t', 'vdj_b', 'vdj_t_gd', 'count', 'antigen_analysis')

FEATURE_TYPES_NAMES = {
            "Gene Expression": "rna",
            "Peaks": "atac",
            "Antibody Capture": "prot",
            "VDJ": "vdj",
            "VDJ-T": "vdj_t",
            "VDJ-B": "vdj_b",
            "CRISPR Guide Capture": "gdo",
            "Multiplexing Capture": "hto",
            "Antigen Capture": "antigen",
        }

def cast_to_writeable_dtype(result: pd.DataFrame) -> pd.DataFrame:
    """
    Cast the dataframe to dtypes that can be written by mudata.
    """    
    # dtype inferral workfs better with np.nan
    result = result.replace({pd.NA: np.nan})

    # MuData supports nullable booleans and ints
    # ie. `IntegerArray` and `BooleanArray`
    result = result.convert_dtypes(infer_objects=True,
                                   convert_integer=True,
                                   convert_string=False,
                                   convert_boolean=True,
                                   convert_floating=False)
    
    # Convert leftover 'object' columns to string
    # However, na values are supported, so convert all values except NA's to string
    object_cols = result.select_dtypes(include='object').columns.values
    for obj_col in object_cols:
        result[obj_col] = result[obj_col].where(result[obj_col].isna(), result[obj_col].astype(str)).astype('category')
    return result

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
    for file_glob in ('*/metrics_summary.csv', '*/count/feature_reference.csv',
                      '*/count/crispr_analysis/perturbation_efficiencies_by_feature.csv',
                      '*/count/crispr_analysis/perturbation_efficiencies_by_target.csv',
                      '*/antigen_analysis',
                     ):
        found_files = list(per_sample_outs_dir.glob(file_glob))
        if len(found_files) > 1:
            raise ValueError(f"Found more than one file for glob '{file_glob}' file. "
                            "This component currently only supports parsing cellranger multi output for one sample.")
        file_name = Path(file_glob).name.removesuffix('.csv')
        found_input[file_name] = found_files[0] if found_files else None

    return found_input


def proces_perturbation(key_name: str, mudata: mudata.MuData, efficiency_file: Path):
    assert 'gdo' in mudata.mod
    eff_df = pd.read_csv(efficiency_file, index_col="Perturbation", sep=",", decimal=".", quotechar='"')
    mudata.mod['gdo'].uns[key_name] = eff_df
    return mudata

def process_feature_reference(mudata: mudata.MuData, efficiency_file: Path):
    df = pd.read_csv(efficiency_file, index_col="id", sep=",", decimal=".", quotechar='"')
    assert 'feature_type' in df.columns, "Columns 'feature_type' should be present in features_reference file."
    feature_types = df['feature_type']
    missing_features = set(feature_types) - set(FEATURE_TYPES_NAMES)
    if missing_features:
        raise ValueError("Not all feature types present in the features_reference file are supported by this component.\n"
                         f"Missing support for features: {','.join(missing_features)}.")
    for feature_type in feature_types:
        modality = FEATURE_TYPES_NAMES[feature_type]
        subset_df = df.loc[df['feature_type'] == feature_type]
        mudata.mod[modality].uns['feature_reference'] = subset_df
    return mudata

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

    feature_types = defaultdict(modality_name_factory, FEATURE_TYPES_NAMES)
    return mudata.MuData(adata, feature_types_names=feature_types)

def process_metrics_summary(mudata: mudata.MuData, metrics_file: Path):
    def read_percentage(val):
        try:
            if str(val).endswith('%'):
                return float(val.strip('%')) / 100
            else:
                return val
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

def process_antigen_analysis(mudata: mudata.MuData, antigen_analysis_folder_path: Path):
    assert 'antigen' in mudata.mod
    per_barcodes_file = antigen_analysis_folder_path / "per_barcode.csv"
    assert per_barcodes_file.is_file(), "Expected a per_barcode.csv file to be present."
    per_barcodes_df = pd.read_csv(per_barcodes_file, index_col="barcode",
                                  sep=",", decimal=".", quotechar='"')
    is_gex_cell = per_barcodes_df['is_gex_cell']
    assert len(set(is_gex_cell.unique().tolist()) - set([False, True])) == 0, \
        "Expected 'is_gex_cell' column to be boolean. Please report this as a bug."
    barcodes_in_gex = per_barcodes_df[is_gex_cell]
    # All of the barcodes listed in the per_barcode.csv with is_gex_cell set to 'True'
    # must be in the 'rna' (an thus also 'antigen') modality
    assert barcodes_in_gex.index.difference(mudata['rna'].obs_names).empty
    orig_obs_names = mudata['antigen'].obs_names.copy() 
    mudata['antigen'].obs = cast_to_writeable_dtype(pd.concat([mudata['antigen'].obs, barcodes_in_gex],
                                                    axis='columns',
                                                    join='outer',
                                                    verify_integrity=True,
                                                    sort=False))
    assert orig_obs_names.equals(mudata['antigen'].obs_names)
    del orig_obs_names

    # The antigen_specificity_scores.csv file is only present when cellranger 
    # multi was run with a [antigen-specificity] section in config
    specificity_file = antigen_analysis_folder_path / "antigen_specificity_scores.csv"
    if specificity_file.is_file():
        antigen_scores_df = pd.read_csv(specificity_file,
                                        index_col=["barcode", "antigen"], sep=",",
                                        decimal=".", quotechar='"')
        score = antigen_scores_df.unstack()
        assert score.index.difference(mudata['rna'].obs_names).empty
        antigens = score.columns.unique(level='antigen')
        for antigen in antigens:
            score_antigen = score.loc[:, (slice(None), antigen)].droplevel("antigen", axis=1)
            score_antigen = score_antigen.reindex(mudata['rna'].obs_names)
            mudata['antigen'].obsm[f'antigen_specificity_scores_{antigen}'] = cast_to_writeable_dtype(score_antigen)
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
        'metrics_summary': process_metrics_summary,
        'feature_reference': process_feature_reference,
        'perturbation_efficiencies_by_feature': partial(proces_perturbation, 'perturbation_efficiencies_by_feature'),
        'perturbation_efficiencies_by_target': partial(proces_perturbation, 'perturbation_efficiencies_by_target'),
        'antigen_analysis': process_antigen_analysis,
    }
    mudata_file = process_counts(input_data['count'])
    for modality_name, modality_data_path in input_data.items():
        if modality_name == "count" or not modality_data_path:
            continue
        try:
            parser_function = dispatcher[modality_name]
        except KeyError as e:
            raise ValueError("This component does not support the "
                             f"parsing of the '{modality_name}' yet.") from e
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