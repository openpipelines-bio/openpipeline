### VIASH START
par = {}
### VIASH END

from scanpy.external.pp import scanorama_integrate
from mudata import read_h5mu

mdata = read_h5mu(par["input"])

mod_name = par["modality"]
mod = mdata.mod[mod_name]

# Integration.
scanorama_integrate(
    mod,
    key=par["obs_batch"],
    basis=par["obsm_input"],
    adjusted_basis=par["obsm_output"],
    knn=par["knn"],
    alpha=par["alpha"],
    sigma=par["sigma"],
    approx=par["approx"],
    batch_size=par["batch_size"],
)

mdata.write_h5mu(par["output"], compression=par["output_compression"])
