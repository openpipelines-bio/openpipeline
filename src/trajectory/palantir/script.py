import sys
import numpy as np
import mudata as mu
import palantir

## VIASH START
par = {
    "input": "resources_test/beyond_test_data/atlas.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca_integrated",
    "start_cell": None,
    "start_cell_cluster": None,
    "start_cell_obs_key": "subpopulation",
    "num_waypoints": 500,
    "n_components": 10,
    "knn": 30,
    "scale_components": True,
    "terminal_states": None,
    "terminal_states_obs_key": None,
    "seed": 42,
    "output": "output.h5mu",
    "obs_pseudotime": "palantir_pseudotime",
    "obs_entropy": "palantir_entropy",
    "obsm_fate_probabilities": "palantir_fate_probabilities",
    "uns_waypoints": "palantir_waypoints",
    "output_compression": None,
}
meta = {
    "cpus": 4,
    "resources_dir": "src/trajectory/palantir/",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _resolve_start_cell(adata, par):
    """Return a single cell barcode to use as the trajectory root."""
    if par["start_cell_cluster"] is not None:
        obs_key = par["start_cell_obs_key"]
        cluster = par["start_cell_cluster"]
        if obs_key not in adata.obs.columns:
            raise ValueError(
                f"--start_cell_obs_key '{obs_key}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        if not (adata.obs[obs_key] == cluster).any():
            raise ValueError(
                f"No cells found for cluster '{cluster}' in obs['{obs_key}']."
            )
        # early_cell(ad, celltype, celltype_column) - requires DM_EigenVectors_multiscaled
        start_cell = palantir.utils.early_cell(adata, cluster, celltype_column=obs_key)
        logger.info(
            "Auto-selected start cell '%s' from cluster '%s' (obs['%s'])",
            start_cell,
            cluster,
            obs_key,
        )
        return start_cell
    elif par["start_cell"] is not None:
        if par["start_cell"] not in adata.obs_names:
            raise ValueError(
                f"--start_cell '{par['start_cell']}' not found in obs_names."
            )
        return par["start_cell"]
    else:
        raise ValueError(
            "Either --start_cell or --start_cell_cluster must be provided."
        )


def _resolve_terminal_states(adata, par):
    """Return a list of terminal-state barcodes, or None for auto-detection."""
    if par["terminal_states"]:
        for cell in par["terminal_states"]:
            if cell not in adata.obs_names:
                raise ValueError(
                    f"--terminal_states barcode '{cell}' not found in obs_names."
                )
        return par["terminal_states"]
    if par["terminal_states_obs_key"]:
        obs_key = par["terminal_states_obs_key"]
        if obs_key not in adata.obs.columns:
            raise ValueError(
                f"--terminal_states_obs_key '{obs_key}' not found in .obs."
            )
        terminal_cells = []
        for label in adata.obs[obs_key].unique():
            mask = adata.obs[obs_key] == label
            sub = adata[mask]
            if "connectivities" in sub.obsp:
                degrees = np.asarray(sub.obsp["connectivities"].sum(axis=1)).ravel()
                best_idx = int(np.argmax(degrees))
            else:
                best_idx = 0
            terminal_cells.append(sub.obs_names[best_idx])
        logger.info(
            "Selected %d terminal states from obs['%s']: %s",
            len(terminal_cells),
            obs_key,
            terminal_cells,
        )
        return terminal_cells
    return None


def main():
    np.random.seed(par["seed"])

    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    obsm_key = par["obsm_input"]
    if obsm_key not in adata.obsm:
        raise ValueError(
            f"--obsm_input '{obsm_key}' not found in .obsm. "
            f"Available keys: {list(adata.obsm.keys())}"
        )
    logger.info(
        "Using embedding '%s' (%d cells, %d dims)",
        obsm_key,
        adata.obsm[obsm_key].shape[0],
        adata.obsm[obsm_key].shape[1],
    )

    # -- 1. Diffusion maps ----------------------------------------------------
    logger.info(
        "Computing diffusion maps (pca_key=%s, n_components=%d, knn=%d)",
        obsm_key,
        par["n_components"],
        par["knn"],
    )
    palantir.utils.run_diffusion_maps(
        adata,
        pca_key=obsm_key,
        n_components=par["n_components"],
        knn=par["knn"],
        seed=par["seed"],
    )

    # -- 2. Multiscale diffusion space (required for early_cell + run_palantir) --
    logger.info("Determining multiscale diffusion space")
    palantir.utils.determine_multiscale_space(adata)

    # -- 3. Resolve start cell ------------------------------------------------
    start_cell = _resolve_start_cell(adata, par)
    logger.info("Start cell: %s", start_cell)

    # -- 4. Resolve terminal states -------------------------------------------
    terminal_states = _resolve_terminal_states(adata, par)
    if terminal_states:
        logger.info("Terminal states (%d): %s", len(terminal_states), terminal_states)
    else:
        logger.info("Terminal states: auto-detected by Palantir")

    # -- 5. Run Palantir ------------------------------------------------------
    # run_palantir stores results directly in adata.obs / adata.obsm / adata.uns
    # when given an AnnData input.
    n_jobs = max(1, (meta.get("cpus") or 1))
    logger.info(
        "Running Palantir (num_waypoints=%d, n_jobs=%d)", par["num_waypoints"], n_jobs
    )
    palantir.core.run_palantir(
        adata,
        early_cell=start_cell,
        terminal_states=terminal_states,
        num_waypoints=par["num_waypoints"],
        scale_components=par["scale_components"],
        n_jobs=n_jobs,
        seed=par["seed"],
        pseudo_time_key=par["obs_pseudotime"],
        entropy_key=par["obs_entropy"],
        fate_prob_key=par["obsm_fate_probabilities"],
        waypoints_key=par["uns_waypoints"],
    )

    logger.info(
        "Results stored: obs['%s'], obs['%s'], obsm['%s'], uns['%s']",
        par["obs_pseudotime"],
        par["obs_entropy"],
        par["obsm_fate_probabilities"],
        par["uns_waypoints"],
    )

    # -- 6. Write output ------------------------------------------------------
    logger.info("Writing output to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    main()
