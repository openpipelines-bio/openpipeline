from typing import Any, Dict, List, Tuple
import math
import tempfile
import subprocess
import tarfile
import gzip
import shutil
from pathlib import Path
import yaml
import pandas as pd
from multiprocess import Pool
import gtfparse
import polars as pl

## VIASH START
par = {
    "input_id": ["mysample1", "mysample2"],
    "input_r1": [
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R1_001.fastq.gz",
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R1_001.fastq.gz",
    ],
    "input_r2": [
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R2_001.fastq.gz",
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R2_001.fastq.gz",
    ],
    "reference_index": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref_v2_7_10_a/",
    "reference_gtf": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz",
    "output": "test_output",
}
meta = {
    "functionality_name": "star_and_htseq",
    "cpus": 30,
    "temp_dir": "/tmp",
    "config": "src/mapping/multi_star/.config.vsh.yaml",
}
## VIASH END

########################
### Helper functions ###
########################


def fetch_arguments_info(config: Dict[str, Any]) -> Dict[str, Any]:
    """Fetch arguments from config"""
    arguments = {
        arg["name"].removeprefix("-").removeprefix("-"): arg
        for group in config["argument_groups"]
        for arg in group["arguments"]
    }
    return arguments

def process_par(
    par: Dict[str, Any],
    arguments_info: Dict[str, Any],
    gz_args: List[str],
    temp_dir: Path
) -> Dict[str, Any]:
    """
    Process the Viash par dictionary

    This turns file strings into Path objects and extracting gzipped files if need be.

    Parameters
    ----------
    par: The par dictionary created by Viash
    arguments_info: The arguments info Dictionary created by `fetch_arguments_info`
    gz_args: A list of argument keys which could be gzip files which need to be decompressed.
    temp_dir: A temporary directory in which to ungzip files
    """
    new_par = {}
    for key, value in par.items():
        arg_info = arguments_info[key]
        # turn file arguments into paths
        if value and arg_info["type"] == "file":
            is_multiple = isinstance(value, list)

            if is_multiple:
                value = [Path(val) for val in value]
            else:
                value = Path(value)

            if key in gz_args:
                print(f">> Checking compression of --{key}", flush=True)
                # turn value into list if need be
                if not is_multiple:
                    value = [value]

                # extract
                value = [extract_if_need_be(path, temp_dir) for path in value]

                # unlist if need be
                if not is_multiple:
                    value = value[0]

        new_par[key] = value
    return new_par

def generate_cmd_arguments(par, arguments_info, step_filter=None, flatten=False):
    """
    Generate command-line arguments by fetching the relevant args

    Parameters
    ----------
    par: The par dictionary created by Viash
    arguments_info: The arguments info Dictionary created by `fetch_arguments_info`
    step_filter: If provided,`par` will be filtered to only contain arguments for which
      argument.info.step == step_filter.
    flatten: If `False`, the command for an argument with multiple values will be
      `["--key", "value1", "--key", "value2"]`, otherwise `["--key", "value1", "value2"]`.
    """
    cmd_args = []

    for key, arg in arguments_info.items():
        arg_val = par.get(key)
        # The info key is always present (changed in viash 0.7.4) 
        # in the parsed config (None if not specified in source config)
        info = arg["info"] or {}
        orig_arg = info.get("orig_arg")
        step = info.get("step")
        if arg_val and orig_arg and (not step_filter or step == step_filter):
            if not arg.get("multiple", False):
                arg_val = [arg_val]

            if arg["type"] in ["boolean_true", "boolean_false"]:
                # if argument is a boolean_true or boolean_false, simply add the flag
                arg_val = [orig_arg]
            elif orig_arg.startswith("-"):
                # if the orig arg flag is not a positional,
                # add the flag in front of each element and flatten
                if flatten:
                    arg_val = [str(x) for x in [orig_arg] + arg_val]
                else:
                    arg_val = [str(x) for val in arg_val for x in [orig_arg, val]]

            cmd_args.extend(arg_val)

    return cmd_args

def is_gz_file(path: Path) -> bool:
    """Check whether something is a gzip"""
    with open(path, "rb") as file:
        return file.read(2) == b"\x1f\x8b"

def extract_if_need_be(par_value: Path, temp_dir_path: Path) -> Path:
    """if {par_value} is a Path, extract it to a temp_dir_path and return the resulting path"""
    if par_value.is_file() and tarfile.is_tarfile(par_value):
        # Remove two extensions (if they exist)
        extaction_dir_name = Path(par_value.stem).stem
        unpacked_path = temp_dir_path / extaction_dir_name
        print(f"  Tar detected; extracting {par_value} to {unpacked_path}", flush=True)

        with tarfile.open(par_value, "r") as open_tar:
            members = open_tar.getmembers()
            root_dirs = [
                member
                for member in members
                if member.isdir() and member.name != "." and "/" not in member.name
            ]
            # if there is only one root_dir (and there are files in that directory)
            # strip that directory name from the destination folder
            if len(root_dirs) == 1:
                for mem in members:
                    mem.path = Path(*Path(mem.path).parts[1:])
            members_to_move = [mem for mem in members if mem.path != Path(".")]
            open_tar.extractall(unpacked_path, members=members_to_move)
        return unpacked_path

    elif par_value.is_file() and is_gz_file(par_value):
        # Remove extension (if it exists)
        extaction_file_name = Path(par_value.stem)
        unpacked_path = temp_dir_path / extaction_file_name
        print(f"  Gzip detected; extracting {par_value} to {unpacked_path}", flush=True)

        with gzip.open(par_value, "rb") as f_in:
            with open(unpacked_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return unpacked_path

    else:
        return par_value

def load_star_reference(reference_index: str) -> None:
    """Load star reference index into memory."""
    subprocess.run(
        [
            "STAR",
            "--genomeLoad", "LoadAndExit",
            "--genomeDir", str(reference_index),
        ],
        check=True
    )

def unload_star_reference(reference_index: str) -> None:
    """Remove star reference index from memory."""
    subprocess.run(
        [
            "STAR",
            "--genomeLoad", "Remove",
            "--genomeDir", str(reference_index),
        ],
        check=True
    )

def star_and_htseq(
    group_id: str,
    r1_files: List[Path],
    r2_files: List[Path],
    temp_dir: Path,
    par: Dict[str, Any],
    arguments_info: Dict[str, Any],
    num_threads: int
) -> Tuple[int, str] :
    star_output = par["output"] / "per" / group_id
    temp_dir_group = temp_dir / f"star_tmp_{group_id}"
    unsorted_bam = star_output / "Aligned.out.bam"
    sorted_bam = star_output / "Aligned.sorted.out.bam"
    counts_file = star_output / "htseq-count.txt"
    multiqc_path = star_output / "multiqc_data"

    print(f">> Running STAR for group '{group_id}' with command:", flush=True)
    star_output.mkdir(parents=True, exist_ok=True)
    temp_dir_group.parent.mkdir(parents=True, exist_ok=True)
    run_star(
        r1_files=r1_files,
        r2_files=r2_files,
        output_dir=star_output,
        temp_dir=temp_dir / f"star_tmp_{group_id}",
        par=par,
        arguments_info=arguments_info,
        num_threads=num_threads
    )
    if not unsorted_bam.exists():
        return (1, f"Could not find unsorted bam at '{unsorted_bam}'")

    if par["run_htseq_count"]:
        print(f">> Running samtools sort for group '{group_id}' with command:", flush=True)
        run_samtools_sort(unsorted_bam, sorted_bam)
        if not sorted_bam.exists():
            return (1, f"Could not find sorted bam at '{unsorted_bam}'")

        print(f">> Running htseq-count for group '{group_id}' with command:", flush=True)
        run_htseq_count(sorted_bam, counts_file, par, arguments_info)
        if not counts_file.exists():
            return (1, f"Could not find counts at '{counts_file}'")

    if par["run_multiqc"]:
        run_multiqc(star_output)
        if not multiqc_path.exists():
            return (1, f"Could not find MultiQC output at '{multiqc_path}'")
    
    return (0, "")

def run_star(
    r1_files: List[Path],
    r2_files: List[Path],
    output_dir: Path,
    temp_dir: Path,
    par: Dict[str, Any],
    arguments_info: Dict[str, Any],
    num_threads: int
) -> None:
    """Run star"""
    # process manual arguments
    r1_pasted = [",".join([str(r1) for r1 in r1_files])]
    r2_pasted = [",".join([str(r2) for r2 in r2_files])] if r2_files else []
    manual_par = {
        "--genomeDir": [par["reference_index"]],
        "--genomeLoad": ["LoadAndRemove"],
        "--runThreadN": [str(num_threads)],
        "--runMode": ["alignReads"],
        "--readFilesIn": r1_pasted + r2_pasted,
        # create a tempdir per group
        "--outTmpDir": [temp_dir],
        # make sure there is a trailing /
        "--outFileNamePrefix": [f"{output_dir}/"],
        # fix the outSAMtype to return unsorted BAM files
        "--outSAMtype": ["BAM", "Unsorted"]
    }
    manual_cmd = [str(x)
        for key, values in manual_par.items()
        for x in [key] + values
    ]

    # process all passthrough star arguments
    par_cmd = generate_cmd_arguments(par, arguments_info, "star", flatten=True)

    # combine into one command and turn into strings
    cmd_args = [str(val) for val in ["STAR"] + manual_cmd + par_cmd]

    # run star
    subprocess.run(cmd_args, check=True)

def run_samtools_sort(
    unsorted_bam: Path,
    sorted_bam: Path
) -> None:
    "Run samtools sort"
    cmd_args = [
        "samtools",
        "sort",
        "-o",
        sorted_bam,
        unsorted_bam,
    ]
    subprocess.run(cmd_args, check=True)

def run_htseq_count(
    sorted_bam: Path,
    counts_file: Path,
    par: Dict[str, Any],
    arguments_info: Dict[str, Any]
) -> None:
    """Run HTSeq count"""
    # process manual arguments
    manual_cmd = [
        sorted_bam,
        par["reference_gtf"]
    ]

    # process all passthrough htseq arguments
    par_cmd = generate_cmd_arguments(par, arguments_info, "htseq")

    # combine into one command and turn into strings
    cmd_args = [str(val) for val in ["htseq-count"] + manual_cmd + par_cmd]

    # run htseq
    with open(counts_file, "w", encoding="utf-8") as file:
        subprocess.run(cmd_args, check=True, stdout=file)

def get_feature_info(reference_gtf) -> pd.DataFrame:
    ref = gtfparse.read_gtf(reference_gtf)
    ref_genes = ref.filter((pl.col("feature") == "gene") | (pl.col("source") == "ERCC"))
    return pd.DataFrame(
        {
            "feature_id": pd.Index(ref_genes.get_column("gene_id")),
            "feature_type": "Gene Expression",
            "feature_name": ref_genes.get_column("gene_name").to_pandas()
        }
    )

def run_multiqc(input_dir: Path) -> None:
    cmd_args = ["multiqc", str(input_dir), "--outdir", str(input_dir), "--no-report", "--force"]

    # run multiqc
    subprocess.run(cmd_args, check=True)


########################
###    Main code     ###
########################

def main(par, meta):
    """Main function"""

    # check input arguments
    assert len(par["input_id"]) == len(par["input_r1"]), "--input_r1 should have same length as --input_id"
    if par["input_r2"]:
        assert len(par["input_id"]) == len(par["input_r2"]), "--input_r2 should have same length as --input_id"

    # read config arguments
    with open(meta["config"], "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # fetch all arguments from the config and turn it into a Dict[str, Argument]
    arguments_info = fetch_arguments_info(config)

    # temp_dir = "tmp/"
    with tempfile.TemporaryDirectory(
        prefix=f"{meta['functionality_name']}-",
        dir=meta["temp_dir"],
        ignore_cleanup_errors=True
    ) as temp_dir:
        temp_dir = Path(temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)

        # turn file strings into Paths and decompress gzip if need be
        gz_args = ["input_r1", "input_r2", "reference_index", "reference_gtf"]
        par = process_par(par, arguments_info, gz_args, temp_dir)

        # make sure input_r2 has same length as input_r1
        if not par["input_r2"]:
            par["input_r2"] = [None for _ in par["input_r1"]]

        # group input_files by input_id
        print(">> Group by --input_id", flush=True)
        grouped_inputs = {}
        for group_id, file_r1, file_r2 in zip(par["input_id"], par["input_r1"], par["input_r2"]):
            if group_id not in grouped_inputs:
                grouped_inputs[group_id] = ([], [])
            grouped_inputs[group_id][0].append(file_r1)
            if file_r2:
                grouped_inputs[group_id][1].append(file_r2)

        # create output dir if need be
        par["output"].mkdir(parents=True, exist_ok=True)

        # store features metadata
        feature_info = get_feature_info(str(par["reference_gtf"]))
        with open(par["output"] / "feature_info.tsv", "w", encoding="utf-8") as file:
            feature_info.to_csv(file, sep="\t", index=False)

        # try:
        #     print(">> Loading genome in memory", flush=True)
        #     load_star_reference(par["reference_index"])

        cpus = meta.get("cpus", 1)
        num_items = len(grouped_inputs)
        pool_size = min(cpus, num_items)
        num_threads_per_task = math.ceil(cpus / pool_size)

        with Pool(pool_size) as pool:
            outs = pool.starmap(
                lambda group_id, files: star_and_htseq(
                    group_id=group_id,
                    r1_files=files[0],
                    r2_files=files[1],
                    temp_dir=temp_dir,
                    par=par,
                    arguments_info=arguments_info,
                    num_threads=num_threads_per_task
                ),
                grouped_inputs.items()
            )

        num_errored = 0
        for exit, msg in outs:
            if exit != 0:
                print(f"Error: {msg}")
                num_errored += 1

        pct_succeeded = 1.0 - num_errored / len(outs)
        print("------------------")
        print(f"Success rate: {math.ceil(pct_succeeded * 100)}%")

        assert pct_succeeded >= par["min_success_rate"], f"Success rate should be at least {math.ceil(par['min_success_rate'] * 100)}%"

if __name__ == "__main__":
    main(par, meta)
