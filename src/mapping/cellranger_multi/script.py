from __future__ import annotations

import sys
import re
import subprocess
import tempfile
import pandas as pd
import yaml
from typing import Optional, Any, Union
import tarfile
from pathlib import Path
import shutil
from itertools import chain

## VIASH START
# The following code has been auto-generated by Viash.
par = {
  'output': './cellranger_test_output',
  'input': ['resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz',
            'resources_test/10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz'],
  'library_id': ['5k_human_antiCMV_T_TBNK_connect_GEX_1_subset',
                 '5k_human_antiCMV_T_TBNK_connect_AB_subset',
                 '5k_human_antiCMV_T_TBNK_connect_VDJ_subset'],
  'library_type': ['Gene Expression', 'Antibody Capture', 'VDJ'],
  'gex_input': None,
  'abc_input': None,
  'cgc_input': None,
  'mux_input': None,
  'vdj_input': None,
  'vdj_t_input': None,
  'vdj_t_gd_input': None,
  'vdj_b_input': None,
  'agc_input': None,
  'library_lanes': None,
  'library_subsample': None,
  'gex_expect_cells': None,
  'gex_chemistry': 'auto',
  'gex_secondary_analysis': False,
  'gex_generate_bam': False,
  'gex_include_introns': False,
  'cell_multiplex_sample_id': None,
  'cell_multiplex_oligo_ids': None,
  'cell_multiplex_description': None,
  'dryrun': False
}
meta = {
  'cpus': 10,
  'memory_b': None,
  'memory_kb': None,
  'memory_mb': None,
  'memory_gb': 15,
  'memory_tb': None,
  'memory_pb': None,
  'temp_dir': '/tmp',
  'config': './target/docker/mapping/cellranger_multi/.config.vsh.yaml',
  'resources_dir': './resources_test'
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

# Tested with cellranger 7.0:
# - omitting the lane number is allowed (e.g. `_L001`)
# - lane number should be omitted across all files if omitted in one
# - replacing `.fastq.` for `.fq.` is NOT allowed
# - omitting `.gz` is allowed

fastq_regex = r'^([A-Za-z0-9\-_\.]+)_S(\d+)_(L(\d+)_)?[RI](\d+)_(\d+)\.fastq(\.gz)?$'
# assert re.match(fastq_regex, "5k_human_GEX_1_subset_S1_L001_R1_001.fastq.gz") is not None
# assert re.match(fastq_regex, "5k_human_GEX_1_subset_S1_R1_001.fastq") is not None
# assert re.match(fastq_regex, "5k_human_GEX_1_subset_S1_R1_001.fastq.gz.txt") is None

# Invert some parameters. Keep the original ones in the config for compatibility
inverted_params = {
    "gex_no_secondary_analysis": "gex_secondary_analysis",
}
for inverted_param, param in inverted_params.items():
    par[inverted_param] = not par[param] if par[param] is not None else None
    del par[param]

GEX_CONFIG_KEYS = {
    "gex_reference": "reference",
    "gex_expect_cells": "expect-cells",
    "gex_force_cells": "force-cells",
    "gex_chemistry": "chemistry",
    "gex_no_secondary_analysis": "no-secondary",
    "gex_generate_bam": "create-bam",
    "gex_include_introns": "include-introns",
    "min_assignment_confidence": "min-assignment-confidence",
    "check_library_compatibility": "check-library-compatibility",
    "barcode_sample_assignment": "barcode-sample-assignment",
    "cmo_set": "cmo-set",
    "probe_set": "probe-set",
    "filter_probes": "filter-probes",
    "gex_r1_length": "r1-length",
    "gex_r2_length": "r2-length",
}

FEATURE_CONFIG_KEYS = {
    "feature_reference": "reference",
    "feature_r1_length": "r1-length",
    "feature_r2_length": "r2-length",
    "min_crispr_umi": "min-crispr-umi",
}

VDJ_CONFIG_KEYS = {"vdj_reference": "reference",
                   "vdj_inner_enrichment_primers": "inner-enrichment-primers",
                   "vdj_r1_length": "r1-length",
                   "vdj_r2_length": "r2-length",
                  }


ANTIGEN_SPECIFICITY_CONFIG_KEYS = {
    "control_id": "control_id",
    "mhc_allele": "mhc_allele",
}


REFERENCE_SECTIONS = {
    "gene-expression": (GEX_CONFIG_KEYS, "index"),
    "feature": (FEATURE_CONFIG_KEYS, "index"),
    "vdj": (VDJ_CONFIG_KEYS, "index"),
    "antigen-specificity": (ANTIGEN_SPECIFICITY_CONFIG_KEYS, "columns"),
}

LIBRARY_CONFIG_KEYS = {'library_id': 'fastq_id',
                       'library_type': 'feature_types',
                       'library_subsample': 'subsample_rate',
                       'library_lanes': 'lanes',
                       'library_chemistry': 'chemistry',
                       }


SAMPLE_PARAMS_CONFIG_KEYS = {'sample_ids': 'sample_id',
                             'cell_multiplex_oligo_ids': 'cmo_ids',
                             'sample_description': 'description',
                             'probe_barcode_ids': 'probe_barcode_ids',
                             'sample_expect_cells': 'expect_cells',
                             'sample_force_cells': 'force_cells'}


# These are derived from the dictionaries above
REFERENCES = tuple(reference_param for reference_param, cellranger_param
                   in chain(GEX_CONFIG_KEYS.items(), FEATURE_CONFIG_KEYS.items(), VDJ_CONFIG_KEYS.items())
                   if cellranger_param == "reference")
LIBRARY_PARAMS = tuple(LIBRARY_CONFIG_KEYS.keys())
SAMPLE_PARAMS = tuple(SAMPLE_PARAMS_CONFIG_KEYS.keys())
HELPER_INPUT = {
    'gex_input': 'Gene Expression',
    'abc_input': 'Antibody Capture',
    'cgc_input': 'CRISPR Guide Capture',
    'mux_input': 'Multiplexing Capture',
    'vdj_input': 'VDJ',
    'vdj_t_input': 'VDJ-T',
    'vdj_t_gd_input': 'VDJ-T-GD',
    'vdj_b_input': 'VDJ-B',
    'agc_input': 'Antigen Capture'
}


def infer_library_id_from_path(input_path: str) -> str:
    match = re.match(fastq_regex, input_path)
    assert match is not None, \
        f"File name of '{input_path}' should match regex {fastq_regex}."
    return match.group(1)

def transform_helper_inputs(par: dict[str, Any]) -> dict[str, Any]:
    helper_input = {
        "input": [],
        "library_id": [],
        "library_type": []
    }
    for input_type, library_type in HELPER_INPUT.items():
        if par[input_type]:
            par[input_type] = resolve_input_directories_to_fastq_paths(par[input_type])

            library_ids = [
                infer_library_id_from_path(path.name) for path in par[input_type]
                ]

            library_id_dict = {}
            for fastq, library_id in zip(par[input_type], library_ids):
                library_id_dict.setdefault(library_id, []).append(fastq)

            for library_id, input in library_id_dict.items():
                helper_input["input"] += input
                helper_input["library_id"].append(library_id)
                helper_input["library_type"].append(library_type)

    assert len(helper_input["library_id"]) == len(set(helper_input["library_id"])), "File names passed to feature type-specific inputs must be unique"

    return helper_input

def lengths_gt1(dic: dict[str, Optional[list[Any]]]) -> dict[str, int]:
    return {key: len(li) for key, li in dic.items()
            if li is not None and isinstance(li, (list, tuple, set))}

def strip_margin(text: str) -> str:
    return re.sub('(\n?)[ \t]*\|', '\\1', text)

def subset_dict(dictionary: dict[str, str],
                keys: Union[dict[str, str], list[str]]) -> dict[str, str]:
    if isinstance(keys, (list, tuple)):
        keys = {key: key for key in keys}
    return {dest_key: dictionary[orig_key]
            for orig_key, dest_key in keys.items()
            if dictionary[orig_key] is not None}

def check_subset_dict_equal_length(group_name: str,
                                   dictionary: dict[str, list[str]]) -> None:
    lens = lengths_gt1(dictionary)
    assert len(set(lens.values())) <= 1, f"The number of values passed to {group_name} "\
                                         f"arguments must be 0, 1 or all the same. Offenders: {lens}"

def resolve_input_directories_to_fastq_paths(input_paths: list[str]) -> list[Path]:

    input_paths = [Path(fastq) for fastq in input_paths]
    if len(input_paths) == 1 and input_paths[0].is_dir():
        logger.info("Detected a directory in input paths, "
                    "traversing to see if we can detect any FASTQ files.")
        input_paths = [input_path for input_path in input_paths[0].rglob('*')
                        if re.match(fastq_regex, input_path.name) ]

    # check input fastq files
    for input_path in input_paths:
        assert re.match(fastq_regex, input_path.name) is not None, \
            f"File name of --input '{input_path}' should match regex {fastq_regex}."

    return input_paths

def make_paths_absolute(par: dict[str, Any], config: Path | str):
    with open(config, 'r', encoding="utf-8") as open_viash_config:
        config = yaml.safe_load(open_viash_config)

    arguments = {
        arg["name"].removeprefix("-").removeprefix("-"): arg
        for group in config["functionality"]["argument_groups"]
        for arg in group["arguments"]
    }
    for arg_name, arg in arguments.items():
        if not par.get(arg_name) or arg["type"] != "file":
            continue
        par_value, is_multiple = par[arg_name], arg["multiple"]
        assert is_multiple in (True, False)
        def make_path_absolute(file: str | Path) -> Path:
            logger.info('Making path %s absolute', file)
            return Path(file).resolve()
        
        new_arg = [make_path_absolute(file) for file in par_value] if is_multiple else make_path_absolute(par_value)
        par[arg_name] = new_arg
    return par

def handle_integers_not_set(par: dict[str, Any], viash_config: Path | str) -> str:
    """
    Allow to use `-1` to define a 'not set' value for arguments of `type: integer` with `multiple: true`.
    """
    with open(viash_config, 'r', encoding="utf-8") as open_viash_config:
        config = yaml.safe_load(open_viash_config)

    arguments = {
        arg["name"].removeprefix("-").removeprefix("-"): arg
        for group in config["functionality"]["argument_groups"]
        for arg in group["arguments"]
    }
    for arg_name, arg in arguments.items():
        if not par.get(arg_name) or arg["type"] != "integer":
            continue
        par_value, is_multiple = par[arg_name], arg["multiple"]
        assert is_multiple in (True, False)

        if not is_multiple:
            continue

        def replace_notset_values(integer_value: int) -> int | None:
            return None if integer_value == -1 else integer_value
        
        # Use an extension array to handle "None" values, otherwise int + NA
        # values would be converted to a "float" dtype
        new_arg = pd.array([replace_notset_values(value) for value in par_value], dtype="Int64")
        par[arg_name] = new_arg
    return par

def process_params(par: dict[str, Any], viash_config: Path | str) -> str:

    if par["input"]:
        assert len(par["library_type"]) > 0, "--library_type must be defined when passing input to --input"
        assert len(par["library_id"]) > 0, "--library_id must be defined when passing input to --input"

        # if par["input"] is a directory, look for fastq files
        par["input"] = resolve_input_directories_to_fastq_paths(par["input"])

    # add helper input
    helper_input = transform_helper_inputs(par)
    for key in ["input", "library_id", "library_type"]:
      par[key] = (par[key] if par[key] else []) + helper_input[key]

      assert len(par[key]) > 0, f"Either --{key} or feature type-specific input (e.g. --gex_input, --abc_input, ...) must be defined"

    # check lengths of libraries metadata
    library_dict = subset_dict(par, LIBRARY_PARAMS)
    check_subset_dict_equal_length("Library", library_dict)

    samples_dict = subset_dict(par, SAMPLE_PARAMS)
    check_subset_dict_equal_length("Samples", samples_dict)

    # Allow using -1 to indicate unset integers for arguments
    # that accept multiple integers.
    par = handle_integers_not_set(par, viash_config)

    # use absolute paths
    return make_paths_absolute(par, viash_config)


def generate_csv_category(name: str, args: dict[str, str], orient: str) -> list[str]:
    assert orient in ("index", "columns")
    if not args:
        return []
    title = [ f'[{name}]' ]
    # Which index to include in csv section is based on orientation
    to_csv_args = {"index": (orient=="index"), "header": (orient=="columns")}
    values = [pd.DataFrame.from_dict(args, orient=orient).to_csv(**to_csv_args).strip()]
    return title + values + [""]


def generate_config(par: dict[str, Any], fastq_dir: str) -> str:
    content_list = []
    par["fastqs"] = fastq_dir
    libraries = dict(LIBRARY_CONFIG_KEYS, **{"fastqs": "fastqs"})
    #TODO: use the union (|) operator when python is updated to 3.9
    all_sections = REFERENCE_SECTIONS | {"libraries": (libraries, "columns"), 
                                         "samples": (SAMPLE_PARAMS_CONFIG_KEYS, "columns")}
    for section_name, (section_params, orientation) in all_sections.items():
        reference_pars = subset_dict(par, section_params)
        content_list += generate_csv_category(section_name, reference_pars, orient=orientation)

    return '\n'.join(content_list)

def main(par: dict[str, Any], meta: dict[str, Any]):
    logger.info("  Processing params")
    par = process_params(par, meta['config'])
    logger.info(par)

    # TODO: throw error or else Cell Ranger will
    with tempfile.TemporaryDirectory(prefix="cellranger_multi-",
                                     dir=meta["temp_dir"]) as temp_dir:
        temp_dir_path = Path(temp_dir)
        for reference_par_name in REFERENCES:
            reference = par[reference_par_name]
            logger.info('Looking at %s to check if it needs decompressing', reference)
            if reference and Path(reference).is_file() and tarfile.is_tarfile(reference):
                extaction_dir_name = Path(reference.stem).stem # Remove two extensions (if they exist)
                unpacked_directory = temp_dir_path / extaction_dir_name
                logger.info('Extracting %s to %s', reference, unpacked_directory)

                with tarfile.open(reference, 'r') as open_tar:
                    members = open_tar.getmembers()
                    root_dirs = [member for member in members if member.isdir()
                                 and member.name != '.' and '/' not in member.name]
                    # if there is only one root_dir (and there are files in that directory)
                    # strip that directory name from the destination folder
                    if len(root_dirs) == 1:
                        for mem in members:
                            mem.path = Path(*Path(mem.path).parts[1:])
                    members_to_move = [mem for mem in members if mem.path != Path('.')]
                    open_tar.extractall(unpacked_directory, members=members_to_move)
                par[reference_par_name] = unpacked_directory

        # Creating symlinks of fastq files to tempdir
        input_symlinks_dir = temp_dir_path / "input_symlinks"
        input_symlinks_dir.mkdir()
        for fastq in par['input']:
            destination = input_symlinks_dir / fastq.name
            destination.symlink_to(fastq)

        logger.info("  Creating config file")
        config_content = generate_config(par, input_symlinks_dir)

        logger.info("  Creating Cell Ranger argument")
        temp_id="run"
        proc_pars=["--disable-ui", "--id", temp_id]

        command_line_parameters = {
            "--localcores": meta['cpus'],
            "--localmem": int(meta['memory_gb']) - 2 if meta['memory_gb'] else None,
        }
        for param, param_value in command_line_parameters.items():
            if param_value:
                proc_pars.append(f"{param}={param_value}")

        ## Run pipeline
        if par["dryrun"]:
            par['output'].mkdir(parents=True, exist_ok=True)

            # write config file
            config_file = par['output'] / "config.csv"
            with open(config_file, "w") as f:
                f.write(config_content)
            proc_pars.append(f"--csv={config_file}")

            # display command that would've been used
            cmd = ["cellranger multi"] + proc_pars + ["--csv=config.csv"]
            logger.info("> " + ' '.join(cmd))
        else:
            # write config file to execution directory
            config_file = temp_dir_path / "config.csv"
            with open(config_file, "w") as f:
                f.write(config_content)
            proc_pars.append(f"--csv={config_file}")

            # Already copy config file to output directory
            par['output'].mkdir(parents=True, exist_ok=True)
            with (par['output'] / "config.csv").open('w') as open_config:
                open_config.write(config_content)

            # run process
            cmd = ["cellranger", "multi"] + proc_pars
            logger.info("> " + ' '.join(cmd))
            process_output = subprocess.run(
                cmd,
                cwd=temp_dir,
                check=False,
                capture_output=True
            )

            with (par["output"] / "cellranger_multi.log").open('w') as open_log:
                open_log.write(process_output.stdout.decode('utf-8'))
            try:
                process_output.check_returncode()
            except subprocess.CalledProcessError as e:
                logger.error(e.output.decode('utf-8'))
                print(e.output.decode('utf-8'), flush=True)
                raise e

            # look for output dir file
            tmp_output_dir = temp_dir_path / temp_id / "outs"
            expected_files = {
                Path("multi"): Path.is_dir,
                Path("per_sample_outs"): Path.is_dir,
                Path("config.csv"): Path.is_file,
            }
            for file_path, type_func in expected_files.items():
                output_path = tmp_output_dir / file_path
                if not type_func(output_path):
                    raise ValueError(f"Could not find expected '{output_path}'")

            for output_path in tmp_output_dir.rglob('*'):
                if output_path.name != "config.csv": # Already created
                    shutil.move(str(output_path), par['output'])

if __name__ == "__main__":
    main(par, meta)