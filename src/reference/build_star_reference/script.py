import tempfile
import subprocess
from pathlib import Path
import tarfile
import gzip
import shutil

## VIASH START
par = {
    "genome_fasta": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/fasta/genome.fa",
    "transcriptome_gtf": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz",
    "output": "star_reference_test",
    "genomeSAindexNbases": 7,
}
meta = {"cpus": 8, "temp_dir": "/tmp"}
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
    "genome_fasta": "genomeFastaFiles",
    "output": "genomeDir",
    "transcriptome_gtf": "sjdbGTFfile",
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
par["genomeDir"].mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(prefix="star-", dir=meta["temp_dir"]) as temp_dir:
    # checking for compressed files, ungzip files if need be
    temp_dir_path = Path(temp_dir)
    for par_name in ["genomeFastaFiles", "sjdbGTFfile"]:
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

    print(">> Constructing command", flush=True)
    par["runMode"] = "genomeGenerate"
    par["outTmpDir"] = temp_dir_path / "run"
    if "cpus" in meta and meta["cpus"]:
        par["runThreadN"] = meta["cpus"]

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
