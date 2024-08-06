import os
import re
import subprocess
import tempfile
from typing import Any
import yaml
import shutil

## VIASH START
par = {
    "genome_fasta": ["resources_test/reference_gencodev41_chr1/reference.fa.gz"],
    "gtf": ["resources_test/reference_gencodev41_chr1/reference.gtf.gz"],
    "extra_sequences": None,
    "mitochondrial_contigs": None,
    "filtering_off": False,
    "wta_only_index": False,
    "extra_star_params": "--genomeSAindexNbases 4",
    "reference_archive": "output.tar.gz",
}
meta = {
    "config": "target/nextflow/reference/build_bdrhap2_reference/.config.vsh.yaml",
    "resources_dir": os.path.abspath("src/reference/build_bdrhap_2_reference"),
    "temp_dir": os.getenv("VIASH_TEMP"),
    "memory_mb": None,
    "cpus": None
}
## VIASH END

def clean_arg(argument):
    argument["clean_name"] = re.sub("^-*", "", argument["name"])
    return argument

def read_config(path: str) -> dict[str, Any]:
    with open(path, "r") as f:
        config = yaml.safe_load(f)
    
    config["functionality"]["arguments"] = [
        clean_arg(arg)
        for grp in config["functionality"]["argument_groups"]
        for arg in grp["arguments"]
    ]
    
    return config

def strip_margin(text: str) -> str:
    return re.sub("(\n?)[ \t]*\|", "\\1", text)

def process_params(par: dict[str, Any], config) -> str:
    # check input parameters
    assert par["genome_fasta"], "Pass at least one set of inputs to --genome_fasta."
    assert par["gtf"], "Pass at least one set of inputs to --gtf."
    assert par["reference_archive"].endswith(".tar.gz"), "Output reference_archive must end with .tar.gz."

    # make paths absolute
    for argument in config["functionality"]["arguments"]:
        if par[argument["clean_name"]] and argument["type"] == "file":
            if isinstance(par[argument["clean_name"]], list):
                par[argument["clean_name"]] = [ os.path.abspath(f) for f in par[argument["clean_name"]] ]
            else:
                par[argument["clean_name"]] = os.path.abspath(par[argument["clean_name"]])
    
    return par

def generate_config(par: dict[str, Any], meta, config) -> str:
    content_list = [strip_margin(f"""\
        |#!/usr/bin/env cwl-runner
        |
        |""")]
        
    config_key_value_pairs = []
    for argument in config["functionality"]["arguments"]:
        config_key = (argument.get("info") or {}).get("config_key")
        arg_type = argument["type"]
        par_value = par[argument["clean_name"]]
        if par_value and config_key:
            config_key_value_pairs.append((config_key, arg_type, par_value))

    if meta["cpus"]:
        config_key_value_pairs.append(("Maximum_threads", "integer", meta["cpus"]))

    # print(config_key_value_pairs)

    for config_key, arg_type, par_value in config_key_value_pairs:
        if arg_type == "file":
            str = strip_margin(f"""\
                |{config_key}:
                |""")
            if isinstance(par_value, list):
                for file in par_value:
                    str += strip_margin(f"""\
                        | - class: File
                        |   location: "{file}"
                        |""")
            else:
                str += strip_margin(f"""\
                    |   class: File
                    |   location: "{par_value}"
                    |""")
            content_list.append(str)
        else:
            content_list.append(strip_margin(f"""\
                |{config_key}: {par_value}
                |"""))
            
    ## Write config to file
    return "".join(content_list)

def get_cwl_file(meta: dict[str, Any]) -> str:
    # create cwl file (if need be)
    cwl_file=os.path.join(meta["resources_dir"], "make_rhap_reference_2.2.1_nodocker.cwl")

    return cwl_file

def main(par: dict[str, Any], meta: dict[str, Any]):
    
    config = read_config(meta["config"])
        
    # Preprocess params
    par = process_params(par, config)

    # fetch cwl file
    cwl_file = get_cwl_file(meta)

    # Create output dir if not exists
    outdir = os.path.dirname(par["reference_archive"])
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ## Run pipeline
    with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]) as temp_dir:
        # Create params file
        config_file = os.path.join(temp_dir, "config.yml")
        config_content = generate_config(par, meta, config)
        with open(config_file, "w") as f:
            f.write(config_content)


        cmd = [
            "cwl-runner",
            "--no-container",
            "--preserve-entire-environment",
            "--outdir",
            temp_dir,
            cwl_file,
            config_file
        ]

        env = dict(os.environ)
        env["TMPDIR"] = temp_dir

        print("> " + " ".join(cmd), flush=True)
        _ = subprocess.check_call(
            cmd,
            cwd=os.path.dirname(config_file),
            env=env
        )

        shutil.move(os.path.join(temp_dir, "reference_bd_rhapsody_v2.tar.gz"), par["reference_archive"])

if __name__ == "__main__":
    main(par, meta)
