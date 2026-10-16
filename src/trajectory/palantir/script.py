import sys
import numpy as np
import mudata as mu
import palantir

## VIASH START
par = {
    "input": "resources_test/beyond_test_data/atlas.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca_integrated",
    "start_id": None,
    "start_cluster": None,
    "start_cluster_column": "subpopulation",
    "num_waypoints": 500,
    "n_components": 10,
    "knn": 30,
    "waypoint_knn": 30,
    "scale_components": True,
    "terminal_states": None,
    "terminal_states_column": None,
    "seed": 42,
    "output": "output.h5mu",
    "pseudotime_column": "palantir_pseudotime",
    "entropy_column": "palantir_entropy",
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


def _resolve_start_id(adata, par):
    """Return a single observation identifier to use as the trajectory root."""
    if par["start_id"] is not None and par["start_cluster"] is not None:
        raise ValueError(
            "--start_id and --start_cluster are mutually exclusive; "
            f"got --start_id '{par['start_id']}' and --start_cluster "
            f"'{par['start_cluster']}'. Provide exactly one."
        )
    if par["start_cluster"] is not None:
        obs_key = par["start_cluster_column"]
        cluster = par["start_cluster"]
        if obs_key not in adata.obs.columns:
            raise ValueError(
                f"--start_cluster_column '{obs_key}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        if not (adata.obs[obs_key] == cluster).any():
            raise ValueError(
                f"No observations found for '{cluster}' in .obs column '{obs_key}'."
            )
        # early_cell(ad, celltype, celltype_column) - requires DM_EigenVectors_multiscaled
        start_id = palantir.utils.early_cell(adata, cluster, celltype_column=obs_key)
        logger.info(
            "Auto-selected root '%s' from population '%s' (column '%s')",
            start_id,
            cluster,
            obs_key,
        )
        return start_id
    elif par["start_id"] is not None:
        if par["start_id"] not in adata.obs_names:
            raise ValueError(
                f"--start_id '{par['start_id']}' is not one of the "
                f"{adata.n_obs} observation identifiers."
            )
        return par["start_id"]
    else:
        raise ValueError("Either --start_id or --start_cluster must be provided.")


def _resolve_terminal_states(adata, par):
    """Return a list of terminal-state identifiers, or None for auto-detection."""
    if par["terminal_states"] and par["terminal_states_column"]:
        raise ValueError(
            "--terminal_states and --terminal_states_column are mutually exclusive; "
            "provide at most one."
        )
    if par["terminal_states"]:
        for name in par["terminal_states"]:
            if name not in adata.obs_names:
                raise ValueError(
                    f"--terminal_states '{name}' is not one of the "
                    f"{adata.n_obs} observation identifiers."
                )
        return par["terminal_states"]
    if par["terminal_states_column"]:
        obs_key = par["terminal_states_column"]
        if obs_key not in adata.obs.columns:
            raise ValueError(
                f"--terminal_states_column '{obs_key}' not found in "
                f".obs. Available columns: {list(adata.obs.columns)}"
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
            "Selected %d terminal states from column '%s': %s",
            len(terminal_cells),
            obs_key,
            terminal_cells,
        )
        return terminal_cells
    return None


def _read_input():
    """Return (adata, obsm_key)."""
    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    obsm_key = par["obsm_input"]
    if obsm_key not in adata.obsm:
        raise ValueError(
            f"--obsm_input '{obsm_key}' not found in .obsm. "
            f"Available keys: {list(adata.obsm.keys())}"
        )
    return adata, obsm_key


def main():
    np.random.seed(par["seed"])

    adata, obsm_key = _read_input()

    logger.info(
        "Using embedding '%s' (%d observations, %d dims)",
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

    # -- 3. Resolve trajectory root ------------------------------------------------
    start_id = _resolve_start_id(adata, par)
    logger.info("Trajectory root: %s", start_id)

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
    n_obs = adata.n_obs
    # Palantir caps the waypoints at the number of observations and then builds a
    # --waypoint_knn graph over them, so too few observations fails inside sklearn
    # with a message that does not name either argument.
    n_waypoints = min(par["num_waypoints"], n_obs)
    if n_waypoints <= par["waypoint_knn"]:
        raise ValueError(
            f"--waypoint_knn ({par['waypoint_knn']}) must be smaller than the number "
            f"of waypoints ({n_waypoints} = min(--num_waypoints "
            f"{par['num_waypoints']}, {n_obs} observations))."
        )
    logger.info(
        "Running Palantir (num_waypoints=%d, waypoint_knn=%d, n_jobs=%d)",
        par["num_waypoints"],
        par["waypoint_knn"],
        n_jobs,
    )
    palantir.core.run_palantir(
        adata,
        early_cell=start_id,
        terminal_states=terminal_states,
        knn=par["waypoint_knn"],
        num_waypoints=par["num_waypoints"],
        scale_components=par["scale_components"],
        n_jobs=n_jobs,
        seed=par["seed"],
        pseudo_time_key=par["pseudotime_column"],
        entropy_key=par["entropy_column"],
        fate_prob_key=par["obsm_fate_probabilities"],
        waypoints_key=par["uns_waypoints"],
    )

    # -- 6. Write output ------------------------------------------------------
    logger.info(
        "Results stored: obs['%s'], obs['%s'], obsm['%s'], uns['%s']",
        par["pseudotime_column"],
        par["entropy_column"],
        par["obsm_fate_probabilities"],
        par["uns_waypoints"],
    )
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
