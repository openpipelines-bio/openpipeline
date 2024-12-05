import re
import tempfile
import subprocess
from pathlib import Path
import tarfile
import gzip
import shutil

## VIASH START
par = {
    "input": [
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R1_001.fastq.gz",
        "resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R2_001.fastq.gz",
        # 'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R1_001.fastq.gz',
        # 'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R2_001.fastq.gz'
    ],
    # 'input': [ 'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/' ],
    "reference": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref_v2_7_10_a/",
    "output": "test_output",
}
meta = {"cpus": 8, "temp_dir": "/tmp"}
## VIASH END

########################
### Helper functions ###
########################

# regex for matching R[12] fastq(gz) files
# examples:
# - TSP10_Fat_MAT_SS2_B134171_B115063_Immune_A1_L003_R1.fastq.gz
# - tinygex_S1_L001_I1_001.fastq.gz
fastqgz_regex = r"(.+)_(R\d+)(_\d+)?\.fastq(\.gz)?"


# helper function for cheching whether something is a gzip
def is_gz_file(path: Path) -> bool:
    with open(path, "rb") as file:
        return file.read(2) == b"\x1f\x8b"


# look for fastq files in a directory
def search_fastqs(path: Path) -> list[Path]:
    if path.is_dir():
        print(
            f"Input '{path}' is a directory, traversing to see if we can detect any FASTQ files.",
            flush=True,
        )
        value_paths = [
            file for file in path.iterdir() if re.match(fastqgz_regex, file.name)
        ]
        return value_paths
    else:
        return [path]


# if {par_value} is a Path, extract it to a temp_dir_path and return the resulting path
def extract_if_need_be(par_value: Path, temp_dir_path: Path) -> Path:
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


########################
###    Main code     ###
########################

# rename keys and convert path strings to Path
# note: only list file arguments here. if non-file arguments also need to be renamed,
# the `processPar()` generator needs to be adapted
to_rename = {
    "input": "readFilesIn",
    "reference": "genomeDir",
    "output": "outFileNamePrefix",
}


def process_par(orig_par, to_rename):
    for key, value in orig_par.items():
        # rename the key in par based on the `to_rename` dict
        if key in to_rename.keys():
            new_key = to_rename[key]

            # also turn value into a Path
            if isinstance(value, list):
                new_value = [Path(val) for val in value]
            else:
                new_value = Path(value)
        else:
            new_key = key
            new_value = value
        yield new_key, new_value


par = dict(process_par(par, to_rename))

# create output dir if need be
par["outFileNamePrefix"].mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(
    prefix="star-", dir=meta["temp_dir"], ignore_cleanup_errors=True
) as temp_dir:
    print(">> Check whether input files are directories", flush=True)
    new_read_files_in = []
    for path in par["readFilesIn"]:
        new_read_files_in.extend(search_fastqs(path))
    par["readFilesIn"] = new_read_files_in
    print("", flush=True)

    # checking for compressed files, ungzip files if need be
    temp_dir_path = Path(temp_dir)
    for par_name in ["genomeDir", "readFilesIn"]:
        par_values = par[par_name]
        if par_values:
            # turn value into list
            is_multiple = isinstance(par_values, list)
            if not is_multiple:
                par_values = [par_values]

            # output list
            new_values = []
            for par_value in par_values:
                print(
                    f">> Check compression of --{par_name} with value: {par_value}",
                    flush=True,
                )
                new_value = extract_if_need_be(par_value, temp_dir_path)
                new_values.append(new_value)

            # unlist if need be
            if not is_multiple:
                new_values = new_values[0]

            # replace value
            par[par_name] = new_values
    # end ungzipping
    print("", flush=True)

    print("Grouping R1/R2 input files into pairs", flush=True)
    input_grouped = {}
    for path in par["readFilesIn"]:
        key = re.search(fastqgz_regex, path.name).group(2)
        if key not in input_grouped:
            input_grouped[key] = []
        input_grouped[key].append(str(path))
    par["readFilesIn"] = [",".join(val) for val in input_grouped.values()]
    print("", flush=True)

    print(">> Constructing command", flush=True)
    par["runMode"] = "alignReads"
    par["outTmpDir"] = temp_dir_path / "run"
    if "cpus" in meta and meta["cpus"]:
        par["runThreadN"] = meta["cpus"]
    # make sure there is a trailing /
    par["outFileNamePrefix"] = f"{par['outFileNamePrefix']}/"

    cmd_args = ["STAR"]
    for name, value in par.items():
        if value is not None:
            if isinstance(value, list):
                cmd_args.extend(["--" + name] + [str(x) for x in value])
            else:
                cmd_args.extend(["--" + name, str(value)])
    print("", flush=True)

    print(">> Running STAR with command:", flush=True)
    print("+ " + " ".join([str(x) for x in cmd_args]), flush=True)
    print("", flush=True)

    subprocess.run(cmd_args, check=True)
