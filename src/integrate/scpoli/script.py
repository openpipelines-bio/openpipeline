import muon as mu
import scanpy as sc
from scarches.models.scpoli import scPoli
import scvi

scvi.settings.seed = par["scvi_tools_random_seed"]

# Load the h5mu object (Input file location will be provided by Viash)
mu_data = mu.read(par["input"])

# Extract a specific modality as an AnnData object (Modality name will be provided by Viash)
adata = mu_data[par["modality"]]

source_data = adata.copy()

# Build the scPoli reference model
scpoli_model = scPoli(
    adata=source_data,
    condition_keys=par["condition_keys"],
    cell_type_keys=par["cell_type_keys"],
    hidden_layer_sizes=par["hidden_layer_sizes"],
    embedding_dims=par["embedding_dims"],
    inject_condition=['encoder', 'decoder'],
)

# Train the scPoli reference model
scpoli_model.train(
    n_epochs=par["n_epochs"], 
    pretraining_epochs=par["pretraining_epochs"], 
    alpha_epoch_anneal=par["alpha_epoch_anneal"],  
    eta=par["eta"])



if par["output_model"]:
    scpoli_model.save(par["output_model"], overwrite=True)

# Get latent representation of the input data
data_latent = scpoli_model.get_latent(
    source_data,
    mean=True
)

adata.obsm[par["obsm_output"]] = data_latent

mu_data.write_h5mu(par["output"], compression=par["output_compression"])
