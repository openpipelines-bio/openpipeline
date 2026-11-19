import os
import re
import subprocess
import tempfile
from typing import Any
import yaml
import shutil
import glob

## VIASH START
par = {
    "genome_fasta": ["resources_test/reference_gencodev41_chr1/reference.fa.gz"],
    "gtf": ["resources_test/reference_gencodev41_chr1/reference.gtf.gz"],
    "extra_sequences": None,
    "transcription_factor_motif_pfm": None,
    "mitochondrial_contigs": None,
    "filtering_off": False,
    "disable_readthrough_filtering": False,
    "filter_pars": False,
    "combine_same_gene_names": False,
    "wta_only_index": False,
    "disable_fasta_validation": False,
    "repair_gtf": False,
    "gtf_validation_complete_mode": False,
    "extra_star_params": "--genomeSAindexNbases 4",
    "reference_archive": "output.tar.gz",
}
meta = {
    "config": "target/nextflow/reference/build_bdrhap_reference/.config.vsh.yaml",
    "resources_dir": os.path.abspath("src/reference/build_bdrhap_reference"),
    "temp_dir": os.getenv("VIASH_TEMP"),
    "memory_mb": None,
    "cpus": None,
}
## VIASH END


def clean_arg(argument):
    argument["clean_name"] = re.sub("^-*", "", argument["name"])
    return argument


def read_config(path: str) -> dict[str, Any]:
    with open(path, "r") as f:
        config = yaml.safe_load(f)

    config["arguments"] = [
        clean_arg(arg) for grp in config["argument_groups"] for arg in grp["arguments"]
    ]

    return config


def process_params(par: dict[str, Any], config) -> str:
    # check input parameters
    assert par["genome_fasta"], "Pass at least one set of inputs to --genome_fasta."
    assert par["gtf"], "Pass at least one set of inputs to --gtf."
    assert par["reference_archive"].endswith(".gz"), (
        "Output reference_archive must end with .tar.gz."
    )

    # make paths absolute
    for argument in config["arguments"]:
        if par[argument["clean_name"]] and argument["type"] == "file":
            if isinstance(par[argument["clean_name"]], list):
                par[argument["clean_name"]] = [
                    os.path.abspath(f) for f in par[argument["clean_name"]]
                ]
            else:
                par[argument["clean_name"]] = os.path.abspath(
                    par[argument["clean_name"]]
                )

    return par


def generate_config(par: dict[str, Any], meta, config) -> str:
    job_inputs = {}
    for argument in config["arguments"]:
        config_key = (argument.get("info") or {}).get("config_key")
        value = par[argument["clean_name"]]
        if not (value and config_key):
            continue
        if argument["type"] == "file":
            if isinstance(value, list):
                value = [{"class": "File", "location": file} for file in value]
            else:
                value = {"class": "File", "location": value}
        job_inputs[config_key] = value

    if meta["cpus"]:
        job_inputs["Maximum_threads"] = int(meta["cpus"])

    return yaml.safe_dump(job_inputs, sort_keys=False)


def get_cwl_file(meta: dict[str, Any]) -> str:
    # create cwl file (if need be)
    cwl_file = os.path.join(meta["resources_dir"], "make_rhap_reference_4.0.cwl")

    return os.path.abspath(cwl_file)


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
    with tempfile.TemporaryDirectory(
        prefix="cwl-bd_rhapsody_wta-", dir=meta["temp_dir"]
    ) as temp_dir:
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
            # temporary directories and intermediate outputs of the pipeline steps
            "--tmpdir-prefix",
            os.path.join(temp_dir, ""),
            cwl_file,
            config_file,
        ]

        print("> " + " ".join(cmd), flush=True)
        _ = subprocess.check_call(cmd, cwd=os.path.dirname(config_file))

        archives = glob.glob(os.path.join(temp_dir, "*.tar.gz"))
        assert len(archives) == 1, (
            f"Expected exactly one reference archive, found: {archives}"
        )
        shutil.move(archives[0], par["reference_archive"])


if __name__ == "__main__":
    main(par, meta)
