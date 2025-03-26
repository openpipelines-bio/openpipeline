import tempfile
import subprocess
from pathlib import Path
import tarfile
import gzip
import shutil
import yaml

## VIASH START
par = {
    "input": ["resources_test/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"],
    "reference": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz",
    "output": "test_output",
}
meta = {"cpus": 2, "temp_dir": "/tmp", "config": "src/mapping/htseq/config.vsh.yaml"}
## VIASH END

########################
### Helper functions ###
########################


# helper function for cheching whether something is a gzip
def is_gz_file(path: Path) -> bool:
    with open(path, "rb") as file:
        return file.read(2) == b"\x1f\x8b"


# if {par_value} is a Path, extract it to a temp_dir_path and return the resulting path
def extract_if_need_be(par_value: Path, temp_dir_path: Path) -> Path:
    if par_value.is_file() and tarfile.is_tarfile(par_value):
        # Remove two extensions (if they exist)
        extaction_dir_name = Path(par_value.stem).stem
        unpacked_path = temp_dir_path / extaction_dir_name
        print(
            f"    Tar detected; extracting {par_value} to {unpacked_path}", flush=True
        )

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
        print(
            f"    Gzip detected; extracting {par_value} to {unpacked_path}", flush=True
        )

        with gzip.open(par_value, "rb") as f_in:
            with open(unpacked_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return unpacked_path

    else:
        return par_value


def generate_args(par, config):
    # fetch arguments from config
    arguments = [
        arg for group in config["argument_groups"] for arg in group["arguments"]
    ]

    cmd_args = []

    for arg in arguments:
        arg_val = par.get(arg["name"].removeprefix("--"))
        orig_arg = arg.get("info", {}).get("orig_arg")
        if arg_val and orig_arg:
            if not arg.get("multiple", False):
                arg_val = [arg_val]

            if arg["type"] in ["boolean_true", "boolean_false"]:
                # if argument is a boolean_true or boolean_false, simply add the flag
                arg_val = [orig_arg]
            elif orig_arg.startswith("-"):
                # if the orig arg flag is not a positional,
                # add the flag in front of each element and flatten
                arg_val = [str(x) for val in arg_val for x in [orig_arg, val]]

            cmd_args.extend(arg_val)

    return cmd_args


########################
###    Main code     ###
########################

# read config arguments
config = yaml.safe_load(Path(meta["config"]).read_text())


with tempfile.TemporaryDirectory(prefix="htseq-", dir=meta["temp_dir"]) as temp_dir:
    # checking for compressed files, ungzip files if need be
    temp_dir_path = Path(temp_dir)
    reference = Path(par["reference"])

    print(f">> Check compression of --reference with value: {reference}", flush=True)
    par["reference"] = extract_if_need_be(reference, temp_dir_path)

    print(">> Constructing command", flush=True)
    cmd_args = ["htseq-count"] + generate_args(par, config)

    # manually process cpus parameter
    if "cpus" in meta and meta["cpus"]:
        cmd_args.extend(["--nprocesses", str(meta["cpus"])])

    print(">> Running htseq-count with command:", flush=True)
    print("+ " + " ".join([str(x) for x in cmd_args]), flush=True)

    subprocess.run(cmd_args, check=True)
