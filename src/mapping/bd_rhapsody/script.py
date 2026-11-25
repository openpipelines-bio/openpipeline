import json
import os
import sys
import re
import subprocess
import tempfile
from typing import Any
import yaml
import shutil
import glob
import psutil

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
    "multiple_alignment_handling": None,
    "multiple_alignment_gene_selection": None,
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
    "biotype_counts": None,
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

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


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


def process_params(par: dict[str, Any], config, temp_dir: str) -> str:
    # check input parameters
    assert par["reads"] or par["reads_atac"], (
        "Pass at least one set of inputs to --reads or --reads_atac."
    )

    # output to temp dir if output_dir was not passed
    if not par["output_dir"]:
        par["output_dir"] = os.path.join(temp_dir, "output")

    # checking sample prefix
    if par["run_name"] and re.match("[^A-Za-z0-9]", par["run_name"]):
        logger.warning(
            "--run_name should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'."
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


def generate_config(par: dict[str, Any], config, meta: dict[str, Any]) -> str:
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

    # Limit the number of threads used within each pipeline step
    if meta["cpus"]:
        job_inputs["Maximum_Threads"] = int(meta["cpus"])

    return yaml.safe_dump(job_inputs, sort_keys=False)


def generate_config_file(
    par: dict[str, Any], config: dict[str, Any], meta: dict[str, Any], temp_dir: str
) -> str:
    config_file = os.path.join(temp_dir, "config.yml")
    config_content = generate_config(par, config, meta)
    with open(config_file, "w") as f:
        f.write(config_content)
    return config_file


def cap_resource(value: float | str, budget: int) -> float | str:
    """Cap a ResourceRequirement value, which is either a number or a CWL expression."""
    if isinstance(value, (int, float)):
        return min(value, budget)
    expression = value.strip() if isinstance(value, str) else ""
    if expression.startswith("${") and expression.endswith("}"):
        body = expression[2:-1]
        return f"${{ return Math.min((function() {{ {body} }})(), {budget}); }}"
    return budget


def cap_resource_requirements(node: Any, memory: int, cpus: int):
    """Recursively cap all ResourceRequirements in a CWL document."""
    if isinstance(node, dict):
        if node.get("class") == "ResourceRequirement":
            if "ramMin" in node:
                node["ramMin"] = cap_resource(node["ramMin"], memory)
            if "coresMin" in node:
                node["coresMin"] = cap_resource(node["coresMin"], cpus)
        for child in node.values():
            cap_resource_requirements(child, memory, cpus)
    elif isinstance(node, list):
        for child in node:
            cap_resource_requirements(child, memory, cpus)


# Memory (in MB) kept for cwltool, the expressions it evaluates and the OS
MEMORY_RESERVE_MB = 2000


def get_resource_budget(meta: dict[str, Any]) -> tuple[int, int]:
    """The memory (in MB) and cpus available to the pipeline steps.

    When no memory or cpus are specified, fall back to the resources available
    to this process (see get_host_memory_mb and get_host_cpus).
    """
    total_memory = int(meta["memory_mb"]) if meta["memory_mb"] else get_host_memory_mb()
    memory = total_memory - MEMORY_RESERVE_MB
    if memory <= 0:
        raise ValueError(
            f"Only {total_memory} MB of memory is available, while {MEMORY_RESERVE_MB} MB "
            "is kept for cwltool and the OS. Please provide more memory."
        )
    cpus = int(meta["cpus"]) if meta["cpus"] else get_host_cpus()
    return memory, cpus


def read_cgroup_file(*paths: str) -> list[str] | None:
    """The whitespace-separated fields of the first existing cgroup file."""
    for path in paths:
        try:
            with open(path) as file:
                return file.read().split()
        except OSError:
            continue
    return None


# cgroup v1 reports 'no limit' as a huge number (e.g. 9223372036854771712)
CGROUP_V1_NO_LIMIT = 2**62


def get_host_memory_mb() -> int:
    """Available memory (in MB).

    When the process has a memory limit (e.g. inside a container), only the headroom
    within that limit is used. The memory of the host is only used without a limit:
    psutil reads /proc/meminfo, which shows the memory of the host inside a container.
    """
    limit = read_cgroup_file(
        "/sys/fs/cgroup/memory.max",  # cgroup v2, 'max' means no limit
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",  # cgroup v1
    )
    if limit and limit[0].isdigit() and int(limit[0]) < CGROUP_V1_NO_LIMIT:
        usage = read_cgroup_file(
            "/sys/fs/cgroup/memory.current",
            "/sys/fs/cgroup/memory/memory.usage_in_bytes",
        )
        used = int(usage[0]) if usage and usage[0].isdigit() else 0
        available = int(limit[0]) - used
    else:
        available = psutil.virtual_memory().available
    return int(available / 2**20)


def get_host_cpus() -> int:
    """Available cpus.

    When the process has a cpu quota (e.g. inside a container), only that quota is
    used. Without a quota, the number of cpus the process is allowed to run on is used.
    """
    quota = read_cgroup_file("/sys/fs/cgroup/cpu.max")  # cgroup v2: '<quota> <period>'
    if not quota:
        v1_quota = read_cgroup_file("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        v1_period = read_cgroup_file("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        if v1_quota and v1_period:
            quota = v1_quota + v1_period
    # 'max' (v2) or -1 (v1) means no quota
    if quota and len(quota) == 2 and quota[0].isdigit() and quota[1].isdigit():
        return max(1, int(quota[0]) // int(quota[1]))
    # the cpus this process is allowed to run on
    return len(os.sched_getaffinity(0))


def generate_cwl_file(meta: dict[str, Any], dir: str, memory: int, cpus: int) -> str:
    # create a cwl file with capped resource requirements
    orig_cwl_file = os.path.abspath(
        os.path.join(meta["resources_dir"], "rhapsody_pipeline_4.0.cwl")
    )

    # Steps requesting more resources than available are never started by cwltool.
    # Cap the requirements of each step to the resources available to this component.
    with open(orig_cwl_file, "r") as file:
        cwl_text = file.read()
    cwl_data = json.loads(cwl_text.removeprefix("#!/usr/bin/env cwl-runner\n"))

    cap_resource_requirements(cwl_data, memory, cpus)

    cwl_file = os.path.join(dir, "pipeline.cwl")
    with open(cwl_file, "w") as file:
        json.dump(cwl_data, file, indent=4)

    return cwl_file


# Output files of the BD Rhapsody pipeline, as glob patterns relative to the output directory
OUTPUT_FILES = {
    "output_seurat": "{run_name}_Seurat.rds",
    "output_mudata": "{run_name}.cellismo",
    "metrics_summary": "{run_name}_Metrics_Summary.csv",
    "pipeline_report": "{run_name}_Pipeline_Report.html",
    "rsec_mols_per_cell": "{run_name}_RSEC_MolsPerCell_MEX.zip",
    "dbec_mols_per_cell": "{run_name}_DBEC_MolsPerCell_MEX.zip",
    "rsec_mols_per_cell_unfiltered": "{run_name}_RSEC_MolsPerCell_Unfiltered_MEX.zip",
    # the BAM files can be prefixed, e.g. 'Combined_sample_Bioproduct.bam'
    "bam": "*{run_name}_Bioproduct.bam",
    "bam_index": "*{run_name}_Bioproduct.bam.*",
    "bioproduct_stats": "{run_name}_Bioproduct_Stats.csv",
    "biotype_counts": "{run_name}_biotypeCounts.csv",
    "dimred_tsne": "{run_name}_tSNE_coordinates.csv",
    "dimred_umap": "{run_name}_UMAP_coordinates.csv",
    "immune_cell_classification": "{run_name}_cell_type_experimental.csv",
    "sample_tag_metrics": "{run_name}_Sample_Tag_Metrics.csv",
    "sample_tag_calls": "{run_name}_Sample_Tag_Calls.csv",
    "sample_tag_counts": "{run_name}_Sample*Tag*.zip",
    "sample_tag_counts_unassigned": "{run_name}_Multiplet_and_Undetermined.zip",
    "vdj_metrics": "{run_name}_VDJ_metrics.csv",
    "vdj_per_cell": "{run_name}_VDJ_perCell.csv",
    "vdj_per_cell_uncorrected": "{run_name}_VDJ_perCell_uncorrected.csv.gz",
    "vdj_dominant_contigs": "{run_name}_VDJ_Dominant_Contigs_AIRR.tsv",
    "vdj_unfiltered_contigs": "{run_name}_VDJ_Unfiltered_Contigs_AIRR.tsv.gz",
    "atac_metrics": "{run_name}_ATAC_Metrics.csv",
    "atac_metrics_json": "{run_name}_ATAC_Metrics.json",
    "atac_fragments": "{run_name}_ATAC_Fragments.bed.gz",
    "atac_fragments_index": "{run_name}_ATAC_Fragments.bed.gz.tbi",
    "atac_transposase_sites": "{run_name}_ATAC_Transposase_Sites.bed.gz",
    "atac_transposase_sites_index": "{run_name}_ATAC_Transposase_Sites.bed.gz.tbi",
    "atac_peaks": "{run_name}_ATAC_Peaks.bed.gz",
    "atac_peaks_index": "{run_name}_ATAC_Peaks.bed.gz.tbi",
    "atac_peak_annotation": "{run_name}_peak_annotation.tsv.gz",
    "atac_cell_by_peak": "{run_name}_ATAC_Cell_by_Peak_MEX.zip",
    "atac_cell_by_peak_unfiltered": "{run_name}_ATAC_Cell_by_Peak_Unfiltered_MEX.zip",
    "atac_bam": "*{run_name}_ATAC.bam",
    "atac_bam_index": "*{run_name}_ATAC.bam.*",
    "protein_aggregates_experimental": "{run_name}_Protein_Aggregates_Experimental.csv",
}


def copy_outputs(par: dict[str, Any]):
    for argument, pattern in OUTPUT_FILES.items():
        destination = par[argument]
        if not destination:
            continue
        pattern = pattern.format(run_name=par["run_name"])
        files = sorted(glob.glob(os.path.join(par["output_dir"], pattern)))
        if isinstance(destination, list):  # multiple: true
            # replace the '*' with the file name without run name and extension
            # e.g. 'sample_SampleTag01_hs.zip' -> 'SampleTag01_hs'
            suffix = pattern.rsplit("*", 1)[1]
            for file in files:
                name = (
                    os.path.basename(file)
                    .removeprefix(f"{par['run_name']}_")
                    .removesuffix(suffix)
                )
                shutil.copy(file, destination[0].replace("*", name))
                logger.info("Copied %s to --%s", file, argument)
        elif len(files) == 1:
            shutil.copy(files[0], destination)
            logger.info("Copied %s to --%s", files[0], argument)
        elif not files:
            logger.warning("No output file found for --%s ('%s')", argument, pattern)
        else:
            raise ValueError(f"Multiple output files found for --{argument}: {files}")


def main(par: dict[str, Any], meta: dict[str, Any], temp_dir: str):
    config = read_config(meta["config"])

    # Preprocess params
    par = process_params(par, config, temp_dir)

    ## Process parameters
    # Run cwltool through a wrapper which limits the resources used by the parallel
    # executor to the resources available to this component instead of the host.
    # Determine the budget once: the steps must be capped to the same values that
    # cwltool uses, and the available memory of the host changes over time.
    memory, cpus = get_resource_budget(meta)
    cmd = [
        sys.executable,
        os.path.abspath(os.path.join(meta["resources_dir"], "run_cwltool.py")),
    ]
    cmd += ["--max_ram_mb", str(memory), "--max_cores", str(cpus)]
    cmd += [
        "--no-container",
        "--preserve-entire-environment",
        "--outdir",
        par["output_dir"],
        # temporary directories and intermediate outputs of the pipeline steps
        "--tmpdir-prefix",
        os.path.join(temp_dir, ""),
    ]

    if par["parallel"]:
        cmd.append("--parallel")

    if par["timestamps"]:
        cmd.append("--timestamps")

    # Create cwl file with the resource requirements capped to the budget
    cwl_file = generate_cwl_file(meta, temp_dir, memory, cpus)
    cmd.append(cwl_file)

    # Create params file
    config_file = generate_config_file(par, config, meta, temp_dir)
    cmd.append(config_file)

    # Create output dir if not exists
    if not os.path.exists(par["output_dir"]):
        os.makedirs(par["output_dir"])

    # Run command
    logger.info("> %s", " ".join(cmd))
    with (
        subprocess.Popen(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            cwd=os.path.dirname(config_file),
        ) as process,
        open(os.path.join(par["output_dir"], "bd_rhapsody.log"), "w") as open_log,
    ):
        for line in process.stdout:
            print(line.strip(), flush=True)
            open_log.write(line)
    logger.info("Checking exit code.")
    # Popen objects do not have a check_exitcode method.
    if process.returncode != 0:
        raise RuntimeError("BD Rhapsody pipeline returned a nonzero exitcode!")

    # Copy outputs
    copy_outputs(par)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(
        prefix="cwl-bd_rhapsody-", dir=meta["temp_dir"]
    ) as temp_dir:
        main(par, meta, temp_dir)
