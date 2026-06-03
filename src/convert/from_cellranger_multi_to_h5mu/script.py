from pathlib import Path
import sys
import scanpy
import pandas as pd
import mudata
import anndata
import numpy as np
from scirpy.io import read_10x_vdj
from collections import defaultdict
from functools import partial
import json
import csv
import tempfile


## VIASH START
par = {
    "input": "resources_test/10x_4plex_dtc/processed/10x_4plex_dtc.cellranger_multi.output",
    "output": "*.h5mu",
    "uns_metrics": "metrics_cellranger",
    "output_compression": "gzip",
    "sample_csv": "samples.csv",
    "output_filtered_data": False,
}
meta = {"resources_dir": "./src/utils/"}
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
    "Custom": "custom",
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


# Expected cellranger multi output layout:
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
#           |   +-- sample_filtered_feature_bc_matrix.h5
#           |   +-- antibody_analysis
#           |   +-- crispr_analysis
#           |       +-- perturbation_efficiencies_by_feature.csv
#           |       +-- perturbation_efficiencies_by_target.csv
#           +-- vdj_t (unused)
#           +-- vdj_b (unused)
#           +-- vdj_t_gd (unused)
#           +-- metrics_summary.csv


def validate_input_directory(dir: Path):
    if not dir.is_dir():
        raise ValueError("Specified input is not a directory.")
    if not (dir / "config.csv").is_file():
        logger.warning(
            "Config.csv not found in input directory, this folder might not be a valid cellranger multi output."
        )
    per_sample_outs = dir / "per_sample_outs"
    if not per_sample_outs.is_dir():
        raise ValueError(
            f"Input folder must contain the subfolder {per_sample_outs} please make "
            "sure that the specified input folder is a valid cellranger multi output."
        )


def _list_library_types(dir: Path) -> dict[str, Path]:
    """List the library-type subfolders under multi/ (or the dir itself if there
    is no multi/ layout). Validates that only known names are present.
    """
    multi_dir = dir / "multi"
    library_type_dir = multi_dir if multi_dir.is_dir() else dir
    found = {}
    for library_type in library_type_dir.iterdir():
        if library_type.name in POSSIBLE_LIBRARY_TYPES:
            if multi_dir.is_dir() and not library_type.is_dir():
                logger.warning(
                    "%s is not a directory. Contents of the multi folder "
                    "must be directories to be recognized as valid input data",
                    library_type,
                )
                continue
            found[library_type.name] = library_type
            continue
        if multi_dir.is_dir():
            raise ValueError(
                f"Contents of the 'multi' folder must be found one of the following: {','.join(POSSIBLE_LIBRARY_TYPES)}."
            )
    return found


def _iter_sample_dirs(dir: Path):
    return [s for s in (dir / "per_sample_outs").iterdir() if s.is_dir()]


def detect_count_matrices(dir: Path) -> dict:
    """Detect count-matrix input: the raw count folder and (optionally) the
    multiplexing-analysis folder used to split it, or per-sample
    already-filtered count matrices when --output_filtered_data is set.
    """
    library_types = _list_library_types(dir)
    found = {}
    for name in ("count", "multiplexing_analysis"):
        if name in library_types:
            found[name] = library_types[name]
    found.setdefault("count", dir)

    if par["output_filtered_data"]:
        found["filtered_counts"] = {}
        for samples_dir in _iter_sample_dirs(dir):
            for file_part in (
                "count/sample_filtered_feature_bc_matrix.h5",
                "sample_filtered_feature_bc_matrix.h5",  # Cell Ranger v10
            ):
                found_file = samples_dir / file_part
                if found_file.exists():
                    found["filtered_counts"][samples_dir.name] = found_file
                    break
            else:
                raise ValueError(
                    f"Expected a filtered count matrix under {samples_dir}, "
                    "but none was found. Make sure the input directory is a "
                    "valid cellranger multi output that contains per-sample "
                    "filtered feature-barcode matrices."
                )
    return found


def detect_modality_input(dir: Path, count_dir: Path) -> dict:
    """Detect auxiliary modality inputs that enrich the per-sample mudatas:
    vdj_*, antigen_analysis, feature_reference, metrics_summary, and the
    perturbation_efficiencies_* files.
    """
    library_types = _list_library_types(dir)
    found = {}
    for name in ("vdj_t", "vdj_b", "vdj_t_gd", "antigen_analysis"):
        if name in library_types:
            found[name] = library_types[name]

    feature_reference = count_dir / "feature_reference.csv"
    if feature_reference.is_file():
        found["feature_reference"] = feature_reference

    for samples_dir in _iter_sample_dirs(dir):
        for file_part in (
            "metrics_summary.csv",
            "count/crispr_analysis/perturbation_efficiencies_by_feature.csv",
            "crispr_analysis/perturbation_efficiencies_by_feature.csv",  # Cell Ranger v10
            "count/crispr_analysis/perturbation_efficiencies_by_target.csv",
            "crispr_analysis/perturbation_efficiencies_by_target.csv",  # Cell Ranger v10
            "antigen_analysis",
        ):
            found_file = samples_dir / file_part
            if found_file.exists():
                file_name = found_file.name.removesuffix(".csv")
                found.setdefault(file_name, {})[samples_dir.name] = found_file
    return found


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
    mudatas: dict[str, mudata.MuData], efficiency_file: dict[str, Path]
):
    df_all = pd.read_csv(
        efficiency_file, index_col="id", sep=",", decimal=".", quotechar='"'
    )
    assert "feature_type" in df_all.columns, (
        "Columns 'feature_type' should be present in features_reference file."
    )
    for _, mudata_obj in mudatas.items():
        df = df_all.copy(deep=True)
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


def _modality_name_factory(library_type):
    return ("".join(library_type.replace("-", "_").split())).lower()


def _rename_var_to_gene_ids(adata: anndata.AnnData):
    # set the gene ids as var_names (unique Ensembl IDs); gene symbols, which
    # scanpy.read_10x_h5 uses as var_names by default, are not unique
    adata.var = adata.var.rename_axis("gene_symbol").reset_index().set_index("gene_ids")


def _aggregated_counts_to_per_sample_mudatas(
    adata: anndata.AnnData, multiplexing_info, metrics_files
):
    logger.info("Convert to mudata")
    feature_types = defaultdict(_modality_name_factory, FEATURE_TYPES_NAMES)
    mudata_all_samples = mudata.MuData(adata, feature_types_names=feature_types)
    if multiplexing_info:
        # Get the mapping between the barcode and the sample ID from one of the metrics files
        metrics_file = pd.read_csv(
            list(metrics_files.values())[0],
            decimal=".",
            quotechar='"',
            thousands=",",
            usecols=("Group Name", "Metric Value", "Metric Name"),
            dtype=pd.StringDtype(),
        )
        sample_ids = metrics_file[metrics_file["Metric Name"] == "Sample ID"]
        barcode_sample_mapping = (
            sample_ids.loc[:, ["Group Name", "Metric Value"]]
            .set_index("Group Name")
            .squeeze()
            .to_dict()
        )
        # When probe barcodes are combined with multiplexing barcodes; the probe barcode (BC)
        # is paired with other sample with other sample information using the syntax "{BC}+{antibody}"
        # (or even "{BC}+{antibody}+{crispr guide}"). The cells_per_tag
        # file only encodes the BCs; so we need to strip the other information.
        barcode_sample_mapping = {
            key_.partition("+")[0]: value_
            for key_, value_ in barcode_sample_mapping.items()
        }
        return split_samples(
            mudata_all_samples, multiplexing_info, barcode_sample_mapping
        )
    return {"run": mudata_all_samples}


def process_counts_filtered(filtered_counts: dict[str, Path]):
    # Unlike the raw matrix, per-sample filtered matrices are already
    # demultiplexed by cellranger, so each h5 maps 1:1 to an output mudata.
    feature_types = defaultdict(_modality_name_factory, FEATURE_TYPES_NAMES)
    mudatas = {}
    for sample_name, filtered_h5 in filtered_counts.items():
        logger.info("Reading %s.", filtered_h5)
        adata = scanpy.read_10x_h5(filtered_h5, gex_only=False)
        _rename_var_to_gene_ids(adata)
        mudatas[sample_name] = mudata.MuData(adata, feature_types_names=feature_types)
    return mudatas


def process_counts(counts_folder: Path, multiplexing_info, metrics_files):
    counts_matrix_file = counts_folder / "raw_feature_bc_matrix.h5"
    logger.info("Reading %s.", counts_matrix_file)
    adata = scanpy.read_10x_h5(counts_matrix_file, gex_only=False)
    _rename_var_to_gene_ids(adata)
    return _aggregated_counts_to_per_sample_mudatas(
        adata, multiplexing_info, metrics_files
    )


def split_samples(mudata_obj, multiplexing_analysis_folder, barcode_sample_mapping):
    result = {}
    cells_per_tag_file = multiplexing_analysis_folder / "cells_per_tag.json"
    with cells_per_tag_file.open("r") as open_json:
        sample_cell_mapping = json.load(open_json)

    for barcode, indices in sample_cell_mapping.items():
        if indices:
            sample_mudata = mudata_obj[indices]
            sample = barcode_sample_mapping[barcode]
            result[sample] = sample_mudata.copy()
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
        assert per_barcodes_file.is_file(), (
            "Expected a per_barcode.csv file to be present."
        )
        per_barcodes_df = pd.read_csv(
            per_barcodes_file, index_col="barcode", sep=",", decimal=".", quotechar='"'
        )
        is_gex_cell = per_barcodes_df["is_gex_cell"]
        assert len(set(is_gex_cell.unique().tolist()) - set([False, True])) == 0, (
            "Expected 'is_gex_cell' column to be boolean. Please report this as a bug."
        )
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
        vdj_anndata = anndata.AnnData()
        json_for_sample = [
            entry for entry in json_obj if entry["barcode"] in mudata_obj.obs_names
        ]
        if json_for_sample:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".json") as tfile:
                json.dump(json_for_sample, tfile, indent=4)
                tfile.flush()
                vdj_anndata = read_10x_vdj(tfile.name)
        mudata_obj.mod[vdj_type] = vdj_anndata
    return mudatas


def build_per_sample_mudatas(
    count_input: dict, modality_input: dict
) -> dict[str, mudata.MuData]:
    """Construct the initial per-sample mudatas from the detected count input.

    When per-sample filtered matrices are available, each maps 1:1 to an output
    mudata. Otherwise the raw aggregated matrix is split using the
    multiplexing_analysis folder, which requires the metrics_summary files to
    resolve the barcode -> sample mapping.
    """
    if count_input.get("filtered_counts"):
        return process_counts_filtered(count_input["filtered_counts"])
    return process_counts(
        count_input["count"],
        count_input.get("multiplexing_analysis"),
        modality_input.get("metrics_summary"),
    )


def get_modalities(
    mudatas: dict[str, mudata.MuData], modality_input: dict
) -> dict[str, mudata.MuData]:
    dispatcher = {
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
    for modality_name, modality_data_path in modality_input.items():
        if not modality_data_path:
            continue
        try:
            parser_function = dispatcher[modality_name]
        except KeyError as e:
            raise ValueError(
                "This component does not support the "
                f"parsing of the '{modality_name}' yet."
            ) from e
        mudatas = parser_function(mudatas, modality_data_path)
    return mudatas


def main():
    cellranger_multi_dir = Path(par["input"])
    # TODO: remove when issue https://github.com/viash-io/viash/issues/706 is resolved.
    if isinstance(par["output"], (list, set, tuple)):
        assert len(par["output"]) == 1, (
            "A single output file template should have been provided."
        )
        par["output"] = par["output"][0]
    assert par["output"].count("*") == 1, (
        f"Expected exactly one wildcard character (*) in output "
        f"files template ({par['output']}). Found {par['output'].count('*')}"
    )
    validate_input_directory(cellranger_multi_dir)
    count_input = detect_count_matrices(cellranger_multi_dir)
    modality_input = detect_modality_input(
        cellranger_multi_dir, count_dir=count_input["count"]
    )
    mudatas = build_per_sample_mudatas(count_input, modality_input)
    result = get_modalities(mudatas, modality_input)
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
