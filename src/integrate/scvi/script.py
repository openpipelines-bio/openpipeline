import mudata
import scvi


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "batch_key": "batch", #Is this true? Make this an input variable with default!
    "counts_layer": "counts", # Check if this is standard
    "vae_epochs": 500,
    "scanvi_epochs": 200,
    "obsm_output": "X_emb_scvi",
    "output": "foo.h5mu",
}
### VIASH END


def main():
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod['rna']

    # Subset to HVG
    adata_hvg = adata[:,adata.var['highly_variable']].copy() #prob memory heavy

    # Set up the data
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        batch_key=par['batch_key'],
        layer=par['counts_layer']
    )

    # Set up the model
    vae_uns = scvi.model.SCVI(
        adata_sub,
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

    # Parameterization (move up to par definition)
    early_stopping_kwargs = {
        'early_stopping': True,
        'early_stopping_monitor': 'elbo_validation',
        'early_stopping_patience': 10,
        'early_stopping_min_delta': 0.0,
    }
    plan_kwargs = {
        "reduce_lr_on_plateau": True,
        "lr_patience": 8,
        "lr_factor": 0.1,
    }

    # Train the model
    vae_uns.train(
        max_epochs=par['vae_epochs'],
        plan_kwargs=['plan_kwargs'],
        **early_stopping_kwargs,
        check_val_every_n_epoch=1,
        use_gpu=True, #Is this always allowed?
    )
    #Note: train_size=1.0 should give better results, but then can't do early_stopping on validation set

    # Get the latent output
    adata.obsm[par['obsm_output']] = vae_uns.get_latent_representation()

    # Write output
    # TODO: check if mdata automatically gets updated when adata does
    #       Alternatively: store vae_uns output in mdata directly
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main()
