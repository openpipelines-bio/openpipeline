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
    "reads": [
        "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz",
        "resources_test/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz",
        "resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz",
        "resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz",
    ],
    "reads_atac": None,
    "reference_archive": "reference_gencodev41_chr1.tar.gz",
    "targeted_reference": [],
    "abseq_reference": [
        "resources_test/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta"
    ],
    "supplemental_reference": [],
    "cell_calling_data": "mRNA",
    "cell_calling_bioproduct_algorithm": None,
    "cell_calling_atac_algorithm": None,
    "exact_cell_count": 4900,
    "expected_cell_count": None,
    "exclude_intronic_reads": None,
    "sample_tags_version": None,
    "tag_names": [],
    "vdj_version": None,
    "predefined_atac_peaks": None,
    "run_name": "sample",
    "generate_bam": False,
    "alignment_star_params": None,
    "alignment_bwa_mem2_params": None,
    "parallel": True,
    "timestamps": False,
    "dryrun": False,
    "output_dir": "output_large_op",
    "output_seurat": "seurat.rds",
    "output_mudata": "mudata.h5mu",
    "metrics_summary": "metrics_summary.csv",
    "pipeline_report": "pipeline_report.html",
    "rsec_mols_per_cell": None,
    "dbec_mols_per_cell": None,
    "rsec_mols_per_cell_unfiltered": None,
    "bam": None,
    "bam_index": None,
    "bioproduct_stats": None,
    "dimred_tsne": None,
    "dimred_umap": None,
    "immune_cell_classification": None,
    "sample_tag_metrics": None,
    "sample_tag_calls": None,
    "sample_tag_counts": None,
    "sample_tag_counts_unassigned": None,
    "vdj_metrics": None,
    "vdj_per_cell": None,
    "vdj_per_cell_uncorrected": None,
    "vdj_dominant_contigs": None,
    "vdj_unfiltered_contigs": None,
    "atac_metrics": None,
    "atac_metrics_json": None,
    "atac_fragments": None,
    "atac_fragments_index": None,
    "atac_transposase_sites": None,
    "atac_transposase_sites_index": None,
    "atac_peaks": None,
    "atac_peaks_index": None,
    "atac_peak_annotation": None,
    "atac_cell_by_peak": None,
    "atac_cell_by_peak_unfiltered": None,
    "atac_bam": None,
    "atac_bam_index": None,
    "protein_aggregates_experimental": None,
    "long_reads": None,
    "custom_star_params": None,
    "custom_bwa_mem2_params": None,
    "abseq_umi": None,
    "target_analysis": None,
    "vdj_jgene_evalue": None,
    "vdj_vgene_evalue": None,
    "write_filtered_reads": None,
}
meta = {
    "config": "target/nextflow/mapping/bd_rhapsody/.config.vsh.yaml",
    "resources_dir": os.path.abspath("src/mapping/bd_rhapsody"),
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


def strip_margin(text: str) -> str:
    return re.sub("(\n?)[ \t]*\|", "\\1", text)


def process_params(par: dict[str, Any], config, temp_dir: str) -> str:
    # check input parameters
    assert (
        par["reads"] or par["reads_atac"]
    ), "Pass at least one set of inputs to --reads or --reads_atac."

    # output to temp dir if output_dir was not passed
    if not par["output_dir"]:
        par["output_dir"] = os.path.join(temp_dir, "output")

    # checking sample prefix
    if par["run_name"] and re.match("[^A-Za-z0-9]", par["run_name"]):
        print(
            "--run_name should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'.",
            flush=True,
        )
        par["run_name"] = re.sub("[^A-Za-z0-9\\-]", "-", par["run_name"])

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


def generate_config(par: dict[str, Any], config) -> str:
    content_list = [
        strip_margin("""\
        |#!/usr/bin/env cwl-runner
        |
        |cwl:tool: rhapsody
        |""")
    ]

    for argument in config["arguments"]:
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
                content_list.append(
                    strip_margin(f"""\
                    |{config_key}: {par[argument["clean_name"]]}
                    |""")
                )

    ## Write config to file
    return "".join(content_list)


def generate_config_file(
    par: dict[str, Any], config: dict[str, Any], temp_dir: str
) -> str:
    config_file = os.path.join(temp_dir, "config.yml")
    config_content = generate_config(par, config)
    with open(config_file, "w") as f:
        f.write(config_content)
    return config_file


def generate_cwl_file(meta: dict[str, Any], dir: str) -> str:
    # create cwl file (if need be)
    orig_cwl_file = os.path.join(
        meta["resources_dir"], "rhapsody_pipeline_2.2.1_nodocker.cwl"
    )

    # Inject computational requirements into pipeline
    if meta["memory_mb"] or meta["cpus"]:
        cwl_file = os.path.join(dir, "pipeline.cwl")

        # Read in the file
        with open(orig_cwl_file, "r") as file:
            cwl_data = file.read()

        # Inject computational requirements into pipeline
        if meta["memory_mb"]:
            memory = int(meta["memory_mb"]) - 2000  # keep 2gb for OS
            cwl_data = re.sub(
                '"ramMin": [^\n]*[^,](,?)\n', f'"ramMin": {memory}\\1\n', cwl_data
            )
        if meta["cpus"]:
            cwl_data = re.sub(
                '"coresMin": [^\n]*[^,](,?)\n',
                f'"coresMin": {meta["cpus"]}\\1\n',
                cwl_data,
            )

        # Write the file out again
        with open(cwl_file, "w") as file:
            file.write(cwl_data)
    else:
        cwl_file = orig_cwl_file

    return cwl_file


def copy_outputs(par: dict[str, Any], config: dict[str, Any]):
    for arg in config["arguments"]:
        par_value = par[arg["clean_name"]]
        if par_value and arg["type"] == "file" and arg["direction"] == "output":
            # example template: '[sample_name]_(assay)_cell_type_experimental.csv'
            template = (arg.get("info") or {}).get("template")
            if template:
                template_glob = (
                    template.replace("sample", par["run_name"])
                    .replace("assay", "*")
                    .replace("number", "*")
                )
                files = glob.glob(os.path.join(par["output_dir"], template_glob))
                if len(files) == 0 and arg["required"]:
                    raise ValueError(
                        f"Expected output file '{template_glob}' not found."
                    )
                elif len(files) > 1 and not arg["multiple"]:
                    raise ValueError(
                        f"Expected single output file '{template_glob}', but found multiple."
                    )

                if not arg["multiple"]:
                    try:
                        shutil.copy(files[0], par_value)
                        print(f"Copied {files[0]} to {par_value}")
                    except IndexError:
                        print(f"Unable to copy {template_glob} to {par_value}")
                else:
                    # replace '*' in par_value with index
                    for i, file in enumerate(files):
                        shutil.copy(file, par_value.replace("*", str(i)))


def main(par: dict[str, Any], meta: dict[str, Any], temp_dir: str):
    config = read_config(meta["config"])

    # Preprocess params
    par = process_params(par, config, temp_dir)

    ## Process parameters
    cmd = [
        "cwl-runner",
        "--no-container",
        "--preserve-entire-environment",
        "--outdir",
        par["output_dir"],
    ]

    if par["parallel"]:
        cmd.append("--parallel")

    if par["timestamps"]:
        cmd.append("--timestamps")

    # Create cwl file (if need be)
    cwl_file = generate_cwl_file(meta, temp_dir)
    cmd.append(cwl_file)

    # Create params file
    config_file = generate_config_file(par, config, temp_dir)
    cmd.append(config_file)

    # keep environment variables but set TMPDIR to temp_dir
    env = dict(os.environ)
    env["TMPDIR"] = temp_dir

    # Create output dir if not exists
    if not os.path.exists(par["output_dir"]):
        os.makedirs(par["output_dir"])

    # Run command
    print("> " + " ".join(cmd), flush=True)
    _ = subprocess.check_call(cmd, cwd=os.path.dirname(config_file), env=env)

    # Copy outputs
    copy_outputs(par, config)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(
        prefix="cwl-bd_rhapsody-", dir=meta["temp_dir"]
    ) as temp_dir:
        main(par, meta, temp_dir)
