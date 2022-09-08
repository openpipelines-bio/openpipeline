### VIASH START

par = {
}
### VIASH END

from mudata import read_h5mu
import bbknn


h5mu_data = read_h5mu(par["input"])
for modality_name in par["modality"]:
    modality = h5mu_data.mod[modality_name]
    bbknn.bbknn(modality,
                use_rep= par["obm_representation"],
                batch_key = par["obs_batch"],
                neighbors_within_batch=par["n_neighbours_within_batch"], 
                n_pcs=par["n_pcs"], 
                trim=par["n_trim"])

h5mu_data.write(par["output"], compression = "gzip")