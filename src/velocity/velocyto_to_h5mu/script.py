import sys
import anndata as ad
import mudata as mu
import numpy as np

numpy_module = sys.modules["numpy"]
numpy_module.string_ = np.bytes_
numpy_module.unicode_ = np.str_
sys.modules["numpy"] = numpy_module

## VIASH START
par = {
    "input_loom": "resources_test/rna_velocity/velocyto_processed/cellranger_tiny.loom",
    "input_h5mu": "/home/rcannood/workspace/openpipelines-bio/openpipeline/resources_test/cellranger_tiny_fastq/raw_dataset.h5mu",
    "modality": "rna_velocity",
    "output": "output.h5mu",
    "layer_spliced": "velo_spliced",
    "layer_unspliced": "velo_unspliced",
    "layer_ambiguous": "velo_ambiguous",
}
## VIASH END

print("Parameters:", par, flush=True)

print("Reading AnnData from loom", flush=True)
adata_in = ad.read_loom(par["input_loom"])
adata_in.var_names = adata_in.var["Accession"]

print("Creating clean AnnData", flush=True)
adata = ad.AnnData(
    obs=adata_in.obs[[]],
    var=adata_in.var[[]],
    layers={
        par["layer_spliced"]: adata_in.layers["spliced"],
        par["layer_unspliced"]: adata_in.layers["unspliced"],
        par["layer_ambiguous"]: adata_in.layers["ambiguous"],
    },
)

if par["input_h5mu"]:
    print("Received input h5mu to read", flush=True)
    mdata = mu.read_h5mu(par["input_h5mu"])

    print(f"Storing AnnData in modality {par['modality']}", flush=True)
    mdata.mod[par["modality"]] = adata
else:
    print("Creating h5mu from scratch", flush=True)
    mdata = mu.MuData({par["modality"]: adata})

print("Resulting mudata:", mdata, flush=True)

print("Writing h5mu to file", flush=True)
mdata.write_h5mu(par["output"], compression=par["output_compression"])
