def check_arguments(par):
    # check output .obs predictions
    if not par["output_obs_predictions"]:
        par["output_obs_predictions"] = [
            t + "_pred" for t in par["reference_obs_targets"]
        ]
    assert len(par["output_obs_predictions"]) == len(par["reference_obs_targets"]), (
        f"Number of output_obs_predictions must match number of reference_obs_targets\npar: {par}"
    )

    # check output .obs uncertainty
    if not par["output_obs_probability"]:
        par["output_obs_probability"] = [
            t + "_probability" for t in par["reference_obs_targets"]
        ]
    assert len(par["output_obs_probability"]) == len(par["reference_obs_targets"]), (
        f"Number of output_obs_probability must match number of reference_obs_targets\npar: {par}"
    )

    return par


def get_reference_features(adata_reference, par, logger):
    if par["reference_obsm_features"] is None:
        logger.info("Using .X of reference data")
        train_data = adata_reference.X
    else:
        logger.info(f"Using .obsm[{par['reference_obsm_features']}] of reference data")
        train_data = adata_reference.obsm[par["reference_obsm_features"]]

    return train_data


def get_query_features(adata, par, logger):
    if par["input_obsm_features"] is None:
        logger.info("Using .X of query data")
        query_data = adata.X
    else:
        logger.info(f"Using .obsm[{par['input_obsm_features']}] of query data")
        query_data = adata.obsm[par["input_obsm_features"]]

    return query_data
