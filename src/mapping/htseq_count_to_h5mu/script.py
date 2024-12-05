import tempfile
from pathlib import Path
import tarfile
import gzip
import shutil
import pandas as pd
import mudata as md
import anndata as ad
import polars as pl
import numpy as np
import gtfparse

## VIASH START
par = {
    "input_counts": ["resources_test/cellranger_tiny_fastq/htseq_counts.tsv"],
    "input_id": ["", "bar"],
    "reference": "resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz",
    "output": "test_output.h5mu",
}
meta = {"temp_dir": "/tmp"}
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


print("> combine counts data", flush=True)
counts_data = []

for input_id, input_counts in zip(par["input_id"], par["input_counts"]):
    data = pd.read_table(
        input_counts,
        index_col=0,
        names=["gene_ids", input_id],
        dtype={"gene_ids": "U", input_id: "i"},
    ).transpose()
    counts_data.append(data)

# combine all counts
counts_and_qc = pd.concat(counts_data, axis=0)

print("> split qc", flush=True)
idx = counts_and_qc.columns.str.startswith("_")
qc = counts_and_qc.loc[:, idx]
qc.columns = qc.columns.str.replace("^__", "", regex=True)
counts = counts_and_qc.loc[:, ~idx]

print("> construct var", flush=True)
with tempfile.TemporaryDirectory(prefix="htseq-", dir=meta["temp_dir"]) as temp_dir:
    # checking for compressed files, ungzip files if need be
    temp_dir_path = Path(temp_dir)
    reference = Path(par["reference"])

    print(f">> Check compression of --reference with value: {reference}", flush=True)
    par["reference"] = extract_if_need_be(reference, temp_dir_path)

    # read_gtf only works on str object, not pathlib.Path
    reference = gtfparse.read_gtf(str(par["reference"]))


# This is a polars dataframe, not pandas
reference_genes = reference.filter(
    (pl.col("feature") == "gene") & (pl.col("gene_id").is_in(list(counts.columns)))
).sort("gene_id")

var = pd.DataFrame(
    data={
        "gene_ids": pd.Index(reference_genes.get_column("gene_id")),
        "feature_types": "Gene Expression",
        "gene_symbol": reference_genes.get_column("gene_name").to_pandas(),
    }
).set_index("gene_ids")

print("> construct anndata", flush=True)
adata = ad.AnnData(X=counts, obsm={"qc_htseq": qc}, var=var, dtype=np.int32)

print("> convert to mudata", flush=True)
mdata = md.MuData(adata)

print("> write to file", flush=True)
mdata.write_h5mu(par["output"], compression=par["output_compression"])
