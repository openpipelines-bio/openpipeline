import subprocess
from pathlib import Path
import mudata as md

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

print("> Running command with folder", flush=True)
input = meta["resources_dir"] + "/multi_star/"
output = "test_output.h5mu"

cmd_pars = [
    meta["executable"],
    "--input",
    str(input),
    "--output",
    output,
    "---cpus",
    "2",
    "--output_compression",
    "gzip",
]
subprocess.run(cmd_pars, check=True)

print("> Check if file exists", flush=True)
output_path = Path(output)
assert output_path.is_file()

print("> Check contents", flush=True)
mdata = md.read_h5mu(output)

print(mdata, flush=True)
print("", flush=True)

assert "rna" in mdata.mod, "MuData should contain RNA modality"
assert mdata.n_obs == 1, "MuData should only contain one observation"

adata_rna = mdata.mod["rna"]
assert adata_rna.n_vars > 100, "RNA modality should contain at least 100 variables"
assert "qc_star" in adata_rna.obsm, "RNA modality should contain STAR QC"
assert "qc_htseq" in adata_rna.obsm, "RNA modality should contain htseq QC"

print("> Completed Successfully!", flush=True)
