import sys
import numpy as np
numpy_module = sys.modules['numpy']
numpy_module.string_ = np.bytes_
numpy_module.unicode_ = np.str_
sys.modules['numpy'] = numpy_module

import subprocess
import pathlib
import mudata
import loompy

## VIASH START
meta = {
    'name': './target/native/convert/from_velocyto_to_h5mu/from_velocyto_to_h5mu',
    'resources_dir': './resources_test/'
}
## VIASH END

tiny_fastq = pathlib.Path(meta["resources_dir"]) / "cellranger_tiny_fastq"
input_loom = tiny_fastq / "velocyto.loom"
input_h5mu = tiny_fastq / "raw_dataset.h5mu"
output = pathlib.Path("output.h5mu")

print(f"Running {meta['name']}", flush=True)
subprocess.run(
    args=[
        meta["executable"],
        "--input_loom",
        input_loom,
        "--input_h5mu",
        input_h5mu,
        "--output",
        output,
        "--output_compression", "gzip"
    ],
    check=True
)

print("Checking whether output exists", flush=True)
assert output.is_file()

print("Reading output file", flush=True)
output_data = mudata.read_h5mu(output)

print("Checking contents", flush=True)
assert list(output_data.mod.keys()) == ["rna", "rna_velocity"]

with loompy.connect(input_loom) as ds:
    mshape = output_data.mod['rna_velocity'].shape[::-1]
    lshape = ds.shape
    assert mshape == lshape, \
        f"Expected mudata shape {mshape} to be the same the loom shape {lshape}"