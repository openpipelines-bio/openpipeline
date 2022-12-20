import mudata
import scvi
from torch.cuda import is_available as cuda_is_available
try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": None,
    "obs_batch": "sample_id",
    "var_input": None,
    "output": "foo.h5mu",
    "obsm_output": "X_scvi_integrated",
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "reduce_lr_on_plateau": True,
    "lr_factor": 0.6,
    "lr_patience": 30,
    "epochs": 500}
### VIASH END


def main():
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod[par['modality']]

    if par['var_input']:
        # Subset to HVG
        adata = adata[:,adata.var['var_input']].copy()

    # Set up the data
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=par['obs_batch'],
        layer=par['input_layer']
    )

    # Set up the model
    vae_uns = scvi.model.SCVI(
        adata,
        n_hidden=128, #this is the default
        n_latent=30,
        n_layers=2,
        dropout_rate=0.1, #this is the default
        dispersion='gene', #this is the default
        gene_likelihood='nb',
        use_layer_norm='both',
        use_batch_norm="none",
        encode_covariates=True, #Parameterization for better scArches performance -> maybe don't use this always?
        deeply_inject_covariates=False, #Parameterization for better scArches performance -> maybe don't use this always?
        use_observed_lib_size=False, #When size_factors are not passed
    )

    plan_kwargs = {
        "reduce_lr_on_plateau": par['reduce_lr_on_plateau'],
        "lr_patience": par['lr_patience'],
        "lr_factor": par['lr_factor'],
    }


    # Train the model
    vae_uns.train(
        max_epochs = par['max_epochs'],
        early_stopping=par['early_stopping'],
        early_stopping_monitor=par['early_stopping_monitor'],
        early_stopping_patience=par['early_stopping_patience'],
        early_stopping_min_delta=par['early_stopping_min_delta'],
        plan_kwargs=plan_kwargs,
        check_val_every_n_epoch=1,
        use_gpu=(cuda_is_available() or mps_is_available()),
    )
    #Note: train_size=1.0 should give better results, but then can't do early_stopping on validation set

    # Get the latent output
    adata.obsm[par['obsm_output']] = vae_uns.get_latent_representation()

    mdata.mod[par['modality']] = adata
    mdata.write_h5mu(par['output'].strip(), compression="gzip")

if __name__ == "__main__":
    main()
