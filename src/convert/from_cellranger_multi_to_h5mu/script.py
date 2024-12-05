from pathlib import Path
import sys
import scanpy
import pandas as pd
import mudata
import numpy as np
from scirpy.io import read_10x_vdj
from collections import defaultdict
from functools import partial
import json
import csv
import tempfile


## VIASH START
par = {
    "input": "resources_test/10x_5k_beam/processed/10x_5k_beam.cellranger_multi.output",
    "output": "foo.h5mu",
    "uns_metrics": "metrics_cellranger",
    "output_compression": "gzip",
}
meta = {"resources_dir": "."}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

POSSIBLE_LIBRARY_TYPES = (
    "vdj_t",
    "vdj_b",
    "vdj_t_gd",
    "count",
    "antigen_analysis",
    "multiplexing_analysis",
)

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
    result = result.convert_dtypes(
        infer_objects=True,
        convert_integer=True,
        convert_string=False,
        convert_boolean=True,
        convert_floating=False,
    )

    # Convert leftover 'object' columns to string
    # However, na values are supported, so convert all values except NA's to string
    object_cols = result.select_dtypes(include="object").columns.values
    for obj_col in object_cols:
        result[obj_col] = (
            result[obj_col]
            .where(result[obj_col].isna(), result[obj_col].astype(str))
            .astype("category")
        )
    return result


def gather_input_data(dir: Path):
    # /
    # +-- multi
    # |    +-- count (raw output)
    # |    |    +-- feature_reference.csv
    # |    |    +-- raw_feature_bc_matrix.h5
    # |    +-- vdj_t
    # |    |    +-- all_contig_annotations.json
    # |    +-- vdj_b
    # |    |    +-- all_contig_annotations.json
    # |    +-- vdj_t_gd
    # |    |    +-- all_contig_annotations.json
    # |    +--  multiplexing_analysis
    # |         +-- cells_per_tag.json
    # +-- per_sample_outs (filtered outputs)
    #      +-- example_1
    #           +-- antigen_analysis
    #           |   +-- per_barcode.csv
    #           |   +-- antigen_specificity_scores.csv
    #           +-- count
    #           |   +-- antibody_analysis
    #           |   +-- crispr_analysis
    #           |       +-- perturbation_efficiencies_by_feature.csv
    #           |       +-- perturbation_efficiencies_by_target.csv
    #           +-- vdj_t (unused)
    #           +-- vdj_b (unused)
    #           +-- vdj_t_gd (unused)
    #           +-- metrics_summary.csv

    if not dir.is_dir():
        raise ValueError("Specified input is not a directory.")
    folder_contents = list(dir.iterdir())
    config = dir / "config.csv"
    if config not in folder_contents:
        logger.warning(
            "Config.csv not found in input directory, this folder might not be a valid cellranger multi output."
        )

    required_subfolders = [
        dir / subfolder_name for subfolder_name in ("multi", "per_sample_outs")
    ]
    found_input = {key_: {} for key_ in POSSIBLE_LIBRARY_TYPES}
    for required_subfolder in required_subfolders:
        if required_subfolder not in folder_contents:
            raise ValueError(
                f"Input folder must contain the subfolder {required_subfolder} please make "
                "sure that the specified input folder is a valid cellranger multi output."
            )

    multi_dir = dir / "multi"
    for library_type in multi_dir.iterdir():
        if not library_type.is_dir():
            logger.warning(
                "%s is not a directory. Contents of the multi folder "
                "must be directories to be recognized as valid input data",
                library_type,
            )
            continue
        if library_type.name not in POSSIBLE_LIBRARY_TYPES:
            raise ValueError(
                f"Contents of the 'multi' folder must be found one of the following: {','.join(POSSIBLE_LIBRARY_TYPES)}."
            )

        found_input[library_type.name] = library_type

    per_sample_outs_dir = dir / "per_sample_outs"
    samples_dirs = [
        samplepath
        for samplepath in per_sample_outs_dir.iterdir()
        if samplepath.is_dir()
    ]
    for samples_dir in samples_dirs:
        for file_part in (
            "metrics_summary.csv",
            "count/feature_reference.csv",
            "count/crispr_analysis/perturbation_efficiencies_by_feature.csv",
            "count/crispr_analysis/perturbation_efficiencies_by_target.csv",
            "antigen_analysis",
        ):
            found_file = samples_dir / file_part
            if found_file.exists():
                file_name = found_file.name.removesuffix(".csv")
                found_input.setdefault(file_name, {})[samples_dir.name] = found_file

    return found_input


def proces_perturbation(
    key_name: str, mudatas: dict[str, mudata.MuData], efficiency_files: dict[str, Path]
):
    for sample_name, mudata_obj in mudatas.items():
        efficiency_file = efficiency_files[sample_name]
        assert "gdo" in mudata_obj.mod
        eff_df = pd.read_csv(
            efficiency_file,
            index_col="Perturbation",
            sep=",",
            decimal=".",
            quotechar='"',
        )
        mudata_obj.mod["gdo"].uns[key_name] = eff_df
    return mudatas


def process_feature_reference(
    mudatas: dict[str, mudata.MuData], efficiency_files: dict[str, Path]
):
    for sample, mudata_obj in mudatas.items():
        efficiency_file = efficiency_files[sample]
        df = pd.read_csv(
            efficiency_file, index_col="id", sep=",", decimal=".", quotechar='"'
        )
        assert (
            "feature_type" in df.columns
        ), "Columns 'feature_type' should be present in features_reference file."
        feature_types = df["feature_type"]
        missing_features = set(feature_types) - set(FEATURE_TYPES_NAMES)
        if missing_features:
            raise ValueError(
                "Not all feature types present in the features_reference file are supported by this component.\n"
                f"Missing support for features: {','.join(missing_features)}."
            )
        for feature_type in feature_types:
            modality = FEATURE_TYPES_NAMES[feature_type]
            subset_df = df.loc[df["feature_type"] == feature_type]
            mudata_obj.mod[modality].uns["feature_reference"] = subset_df
    return mudatas


def process_counts(counts_folder: Path, multiplexing_info, metrics_files):
    counts_matrix_file = counts_folder / "raw_feature_bc_matrix.h5"
    logger.info("Reading %s.", counts_matrix_file)
    adata = scanpy.read_10x_h5(counts_matrix_file, gex_only=False)

    # set the gene ids as var_names
    logger.info("Renaming var columns")
    adata.var = adata.var.rename_axis("gene_symbol").reset_index().set_index("gene_ids")

    # generate output
    logger.info("Convert to mudata")

    def modality_name_factory(library_type):
        return ("".join(library_type.replace("-", "_").split())).lower()

    feature_types = defaultdict(modality_name_factory, FEATURE_TYPES_NAMES)
    mudata_all_samples = mudata.MuData(adata, feature_types_names=feature_types)
    if multiplexing_info:
        # Get the mapping between the barcode and the sample ID from one of the metrics files
        metrics_file = pd.read_csv(
            list(metrics_files.values())[0], decimal=".", quotechar='"', thousands=","
        )
        sample_ids = metrics_file[
            (metrics_file["Metric Name"] == "Sample ID")
            & (metrics_file["Grouped By"] == "Probe barcode ID")
        ]
        barcode_sample_mapping = (
            sample_ids.loc[:, ["Group Name", "Metric Value"]]
            .set_index("Group Name")
            .squeeze()
            .to_dict()
        )
        return split_samples(
            mudata_all_samples, multiplexing_info, barcode_sample_mapping
        )
    return {"run": mudata_all_samples}


def split_samples(mudata_obj, multiplexing_analysis_folder, barcode_sample_mapping):
    result = {}
    cells_per_tag_file = multiplexing_analysis_folder / "cells_per_tag.json"
    with cells_per_tag_file.open("r") as open_json:
        sample_cell_mapping = json.load(open_json)

    for barcode, indices in sample_cell_mapping.items():
        if indices:
            sample_mudata = mudata_obj[indices]
            result[barcode_sample_mapping[barcode]] = sample_mudata.copy()
    return result


def process_metrics_summary(
    mudatas: dict[str, mudata.MuData], metrics_files: dict[str, Path]
):
    def read_percentage(val):
        try:
            if str(val).endswith("%"):
                return float(val.strip("%")) / 100
            else:
                return val
        except (AttributeError, ValueError):
            return val

    for sample, mudata_obj in mudatas.items():
        metrics_file = metrics_files[sample]
        metrics_summary = pd.read_csv(
            metrics_file, decimal=".", quotechar='"', thousands=","
        ).applymap(read_percentage)

        mudata_obj.uns[par["uns_metrics"]] = metrics_summary
        for colname, coldata in metrics_summary.items():
            try:
                new_column = coldata.astype(str, copy=True).astype(
                    {colname: "category"}
                )
                metrics_summary[colname] = new_column
            except (ValueError, TypeError):
                logger.warning(f"Could not store column {colname} from metrics.")
                pass
    return mudatas


def process_antigen_analysis(
    mudatas: dict[str, mudata.MuData], antigen_analysis_folder_paths: dict[str, Path]
):
    for sample_id, mudata_obj in mudatas.items():
        antigen_analysis_folder_path = antigen_analysis_folder_paths[sample_id]
        assert "antigen" in mudata_obj.mod
        per_barcodes_file = antigen_analysis_folder_path / "per_barcode.csv"
        assert (
            per_barcodes_file.is_file()
        ), "Expected a per_barcode.csv file to be present."
        per_barcodes_df = pd.read_csv(
            per_barcodes_file, index_col="barcode", sep=",", decimal=".", quotechar='"'
        )
        is_gex_cell = per_barcodes_df["is_gex_cell"]
        assert (
            len(set(is_gex_cell.unique().tolist()) - set([False, True])) == 0
        ), "Expected 'is_gex_cell' column to be boolean. Please report this as a bug."
        barcodes_in_gex = per_barcodes_df[is_gex_cell]
        # All of the barcodes listed in the per_barcode.csv with is_gex_cell set to 'True'
        # must be in the 'rna' (an thus also 'antigen') modality
        assert barcodes_in_gex.index.difference(mudata_obj["rna"].obs_names).empty
        orig_obs_names = mudata_obj["antigen"].obs_names.copy()
        mudata_obj["antigen"].obs = cast_to_writeable_dtype(
            pd.concat(
                [mudata_obj["antigen"].obs, barcodes_in_gex],
                axis="columns",
                join="outer",
                verify_integrity=True,
                sort=False,
            )
        )
        assert orig_obs_names.equals(mudata_obj["antigen"].obs_names)
        del orig_obs_names

        # The antigen_specificity_scores.csv file is only present when cellranger
        # multi was run with a [antigen-specificity] section in config
        specificity_file = (
            antigen_analysis_folder_path / "antigen_specificity_scores.csv"
        )
        if specificity_file.is_file():
            antigen_scores_df = pd.read_csv(
                specificity_file,
                index_col=["barcode", "antigen"],
                sep=",",
                decimal=".",
                quotechar='"',
            )
            score = antigen_scores_df.unstack()
            assert score.index.difference(mudata_obj["rna"].obs_names).empty
            antigens = score.columns.unique(level="antigen")
            for antigen in antigens:
                score_antigen = score.loc[:, (slice(None), antigen)].droplevel(
                    "antigen", axis=1
                )
                score_antigen = score_antigen.reindex(mudata_obj["rna"].obs_names)
                mudata_obj["antigen"].obsm[f"antigen_specificity_scores_{antigen}"] = (
                    cast_to_writeable_dtype(score_antigen)
                )
    return mudatas


def process_vdj(mudatas: dict[str, mudata.MuData], vdj_folder_path: str):
    # https://scverse.org/scirpy/latest/generated/scirpy.io.read_10x_vdj.html#scirpy-io-read-10x-vdj
    # According to docs, using the json is preferred as this file includes intron info.
    all_config_json_file = vdj_folder_path / "all_contig_annotations.json"
    vdj_type = vdj_folder_path.name
    with all_config_json_file.open("r") as open_json:
        json_obj = json.load(open_json)
    for _, mudata_obj in mudatas.items():
        json_for_sample = [
            entry for entry in json_obj if entry["barcode"] in mudata_obj.obs_names
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json") as tfile:
            json.dump(json_for_sample, tfile, indent=4)
            tfile.flush()
            vdj_anndata = read_10x_vdj(tfile.name)
        mudata_obj.mod[vdj_type] = vdj_anndata
    return mudatas


def get_modalities(input_data):
    dispatcher = {
        "multiplexing_analysis": split_samples,
        "vdj_t": process_vdj,
        "vdj_b": process_vdj,
        "vdj_t_gd": process_vdj,
        "metrics_summary": process_metrics_summary,
        "feature_reference": process_feature_reference,
        "perturbation_efficiencies_by_feature": partial(
            proces_perturbation, "perturbation_efficiencies_by_feature"
        ),
        "perturbation_efficiencies_by_target": partial(
            proces_perturbation, "perturbation_efficiencies_by_target"
        ),
        "antigen_analysis": process_antigen_analysis,
    }
    mudata_per_sample = process_counts(
        input_data["count"],
        input_data["multiplexing_analysis"],
        input_data["metrics_summary"],
    )
    for modality_name, modality_data_path in input_data.items():
        if (
            modality_name in ("count", "multiplexing_analysis")
            or not modality_data_path
        ):
            continue
        try:
            parser_function = dispatcher[modality_name]
        except KeyError as e:
            raise ValueError(
                "This component does not support the "
                f"parsing of the '{modality_name}' yet."
            ) from e
        mudata_per_sample = parser_function(mudata_per_sample, modality_data_path)
    return mudata_per_sample


def main():
    cellranger_multi_dir = Path(par["input"])
    # TODO: remove when issue https://github.com/viash-io/viash/issues/706 is resolved.
    if isinstance(par["output"], (list, set, tuple)):
        assert (
            len(par["output"]) == 1
        ), "A single output file template should have been provided."
        par["output"] = par["output"][0]
    assert par["output"].count("*") == 1, (
        f"Expected exactly one wildcard character (*) in output "
        f"files template ({par['output']}). Found {par['output'].count('*')}"
    )
    input_data = gather_input_data(cellranger_multi_dir)
    result = get_modalities(input_data)
    output_files = {
        par["output"].replace("*", sample_name) for sample_name in result.keys()
    }
    assert len(output_files) == len(result.keys()), (
        "Replacing the wildcard in the output files "
        "template did not produce unique file paths."
    )
    logger.info(
        "Writing output for samples: '%s' to '%s'",
        "".join(result.keys()),
        par["output"],
    )
    with Path(par["sample_csv"]).open("w", newline="") as open_csv:
        csvwriter = csv.DictWriter(open_csv, fieldnames=["sample_name", "file"])
        csvwriter.writeheader()
        for sample_name, mudata_obj in result.items():
            output_file = Path(par["output"].replace("*", sample_name))
            mudata_obj.write_h5mu(output_file, compression=par["output_compression"])
            csvwriter.writerow({"sample_name": sample_name, "file": output_file.name})


if __name__ == "__main__":
    main()
