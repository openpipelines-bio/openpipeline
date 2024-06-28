import os
import re
import subprocess
import tempfile
from typing import Any
import yaml

## VIASH START
par = {
    'reads': [
        'resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz', 
        'resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz'
    ],
    'reads_atac': None,
    'reference_archive': "resources_test/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz",
    'targeted_reference': [],
    'abseq_reference': [],
    'supplemental_reference': [],
    'output': 'output_dir',
    'cell_calling_data': None,
    'cell_calling_bioproduct_algorithm': None,
    'cell_calling_atac_algorithm': None,
    'exact_cell_count': None,
    'expected_cell_count': None,
    'exclude_intronic_reads': None,
    'sample_tags_version': None,
    'tag_names': [],
    'vdj_version': None,
    'predefined_atac_peaks': None,
    'run_name': "sample",
    'generate_bam': None,
    'alignment_star_params': None,
    'alignment_bwa_mem2_params': None,
    'parallel': True,
    'timestamps': False,
    'dryrun': False
}
meta = {
    'config': "target/nextflow/mapping/bd_rhapsody_2/.config.vsh.yaml",
    'resources_dir': os.path.abspath('src/mapping/bd_rhapsody_2'),
    'temp_dir': os.getenv("VIASH_TEMP"),
    'memory_mb': None,
    'cpus': None
}
## VIASH END

def clean_arg(argument):
    argument["clean_name"] = re.sub("^-*", "", argument["name"])
    return argument

def read_config(path: str) -> dict[str, Any]:
    with open(path, 'r') as f:
        config = yaml.safe_load(f)
    
    config["functionality"]["arguments"] = [
        clean_arg(arg)
        for grp in config["functionality"]["argument_groups"]
        for arg in grp["arguments"]
    ]
    
    return config

def logger(msg: str):
    print(msg, flush=True)

def strip_margin(text: str) -> str:
    return re.sub('(\n?)[ \t]*\|', '\\1', text)

def process_params(par: dict[str, Any], config) -> str:
    # check input parameters
    assert par["reads"] or par["reads_atac"], "Pass at least one set of inputs to --reads or --reads_atac."

    # checking sample prefix
    if par["run_name"] and re.match("[^A-Za-z0-9]", par["run_name"]):
        logger("--run_name should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'.")
        par["run_name"] = re.sub("[^A-Za-z0-9\\-]", "-", par["run_name"])

    # make paths absolute
    for argument in config["functionality"]["arguments"]:
        if par[argument["clean_name"]] and argument["type"] == "file":
            if isinstance(par[argument["clean_name"]], list):
                par[argument["clean_name"]] = [ os.path.abspath(f) for f in par[argument["clean_name"]] ]
            else:
                par[argument["clean_name"]] = os.path.abspath(par[argument["clean_name"]])
    
    return par

def generate_config(par: dict[str, Any], config) -> str:
    content_list = [strip_margin(f"""\
        |#!/usr/bin/env cwl-runner
        |
        |cwl:tool: rhapsody
        |""")]

    for argument in config["functionality"]["arguments"]:
        arg_info = argument.get("info") or {}
        config_key = arg_info.get("config_key")
        if par[argument["clean_name"]] and config_key:

            if argument["type"] == "file":
                str = strip_margin(f"""\
                    |{config_key}:
                    |""")
                if isinstance(par[argument["clean_name"]], list):
                    for file in par[argument["clean_name"]]:
                        str += strip_margin(f"""\
                            | - class: File
                            |   location: "{file}"
                            |""")
                else:
                    str += strip_margin(f"""\
                        |   class: File
                        |   location: "{par[argument["clean_name"]]}"
                        |""")
                content_list.append(str)
            else:
                content_list.append(strip_margin(f"""\
                    |{config_key}: {par[argument["clean_name"]]}
                    |"""))

    ## Write config to file
    return ''.join(content_list)

def generate_cwl_file(par: dict[str, Any], meta: dict[str, Any]) -> str:
    # create cwl file (if need be)
    orig_cwl_file=os.path.join(meta["resources_dir"], "rhapsody_pipeline_2.2.1_nodocker.cwl")

    # Inject computational requirements into pipeline
    if meta["memory_mb"] or meta["cpus"]:
        cwl_file = os.path.join(par["output"], "pipeline.cwl")

        # Read in the file
        with open(orig_cwl_file, 'r') as file :
            cwl_data = file.read()

        # Inject computational requirements into pipeline
        if meta["memory_mb"]:
            memory = int(meta["memory_mb"]) - 2000 # keep 2gb for OS
            cwl_data = re.sub('"ramMin": [^\n]*[^,](,?)\n', f'"ramMin": {memory}\\1\n', cwl_data)
        if meta["cpus"]:
            cwl_data = re.sub('"coresMin": [^\n]*[^,](,?)\n', f'"coresMin": {meta["cpus"]}\\1\n', cwl_data)

        # Write the file out again
        with open(cwl_file, 'w') as file:
            file.write(cwl_data)
    else:
        cwl_file = orig_cwl_file

    return cwl_file

def main(par: dict[str, Any], meta: dict[str, Any]):
    config = read_config(meta["config"])
        
    # Preprocess params
    par = process_params(par, config)

    # Create output dir if not exists
    if not os.path.exists(par["output"]):
        os.makedirs(par["output"])

    ## Process parameters
    proc_pars = [
        "--no-container",
        "--preserve-entire-environment",
        "--outdir",
        par["output"]
    ]

    if par["parallel"]:
        proc_pars.append("--parallel")

    if par["timestamps"]:
        proc_pars.append("--timestamps")

    # Create params file
    config_file = os.path.join(par["output"], "config.yml")
    config_content = generate_config(par, config)
    with open(config_file, "w") as f:
        f.write(config_content)

    # Create cwl file (if need be)
    cwl_file = generate_cwl_file(par, meta)

    ## Run pipeline
    if not par["dryrun"]:
        with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]) as temp_dir:
            cmd = ["cwl-runner"] + proc_pars + [cwl_file, config_file]

            env = dict(os.environ)
            env["TMPDIR"] = temp_dir

            logger("> " + ' '.join(cmd))
            _ = subprocess.check_call(
                cmd,
                cwd=os.path.dirname(config_file),
                env=env
            )

    # maybe this won't be necessary anymore
    # # extracting feature ids from references
    # # extract info from reference files (while they still exist)
    # feature_df = extract_feature_types(par)
    # feature_types_file = os.path.join(par["output"], "feature_types.tsv")
    # feature_df.to_csv(feature_types_file, sep="\t", index=False)

    # if not par["dryrun"]:
    #     # look for counts file
    #     if not par["run_name"]:
    #         par["run_name"] = "sample"
    #     counts_filename = par["run_name"] + "_RSEC_MolsPerCell.csv"
        
    #     if par["sample_tags_version"]:
    #         counts_filename = "Combined_" + counts_filename
    #     counts_file = os.path.join(par["output"], counts_filename)
        
    #     if not os.path.exists(counts_file):
    #         raise ValueError(f"Could not find output counts file '{counts_filename}'")

    #     # look for metrics file
    #     metrics_filename = par["run_name"] + "_Metrics_Summary.csv"
    #     metrics_file = os.path.join(par["output"], metrics_filename)
    #     if not os.path.exists(metrics_file):
    #         raise ValueError(f"Could not find output metrics file '{metrics_filename}'")

if __name__ == "__main__":
    main(par, meta)
