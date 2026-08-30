import sys
import scanpy as sc
import mudata as mu

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "obs_groups": "harmony_integration_leiden_1.0",
    "uns_neighbors": "harmonypy_integration_neighbors",
    "use_rna_velocity": False,
    "uns_velocity_graph": "velocity_graph",
    "model": "v1.2",
    "copy": False,
    "output": "output.h5mu",
    "uns_output": "paga",
    "output_compression": None,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading %s, modality %s", par["input"], par["modality"])
data = mu.read_h5ad(par["input"], mod=par["modality"])

if par["obs_groups"] not in data.obs.columns:
    raise ValueError(
        f"Requested to use .obs column {par['obs_groups']} as the grouping "
        f"for PAGA, but the column is not available for modality {par['modality']}."
    )

if par["uns_neighbors"] not in data.uns:
    raise ValueError(
        f"Requested to use .uns key {par['uns_neighbors']} for the neighbors "
        f"settings, but the key is not available for modality {par['modality']}."
    )

# scanpy's RNA-velocity branch always reads the graph from .uns['velocity_graph'],
# so when a different key is requested, alias its contents there for the duration
# of the call and remove the alias again afterwards.
velocity_graph_key = par["uns_velocity_graph"]
velocity_graph_aliased = False
if par["use_rna_velocity"] and velocity_graph_key != "velocity_graph":
    if velocity_graph_key not in data.uns:
        raise ValueError(
            f"Requested to use .uns key {velocity_graph_key} for the RNA velocity "
            f"graph, but the key is not available for modality {par['modality']}."
        )
    data.uns["velocity_graph"] = data.uns[velocity_graph_key]
    velocity_graph_aliased = True

logger.info("Running PAGA")
result = sc.tl.paga(
    data,
    groups=par["obs_groups"],
    use_rna_velocity=par["use_rna_velocity"],
    model=par["model"],
    neighbors_key=par["uns_neighbors"],
    copy=par["copy"],
)
if par["copy"]:
    data = result

if velocity_graph_aliased:
    del data.uns["velocity_graph"]

if par["uns_output"] != "paga":
    data.uns[par["uns_output"]] = data.uns.pop("paga")

logger.info("Writing output to %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    output_file=par["output"],
    h5mu=par["input"],
    modality_name=par["modality"],
    modality_data=data,
    output_compression=par["output_compression"],
)
