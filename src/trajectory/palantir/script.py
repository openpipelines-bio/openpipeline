import sys
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
import palantir

## VIASH START
par = {
    "input": "resources_test/beyond_test_data/atlas.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca_integrated",
    "input_table": None,
    "id_column": None,
    "metadata": None,
    "start_group": None,
    "start_group_cluster": None,
    "start_group_column": "subpopulation",
    "num_waypoints": 500,
    "n_components": 10,
    "knn": 30,
    "waypoint_knn": 30,
    "scale_components": True,
    "terminal_states": None,
    "terminal_states_column": None,
    "seed": 42,
    "output": "output.h5mu",
    "output_table": None,
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
from group_table import read_group_table, write_group_table

logger = setup_logger()

# Key of the embedding inside the in-memory AnnData built in table mode. Internal
# only: table mode never reads or writes an .obsm.
TABLE_EMBEDDING_KEY = "X_embedding"


def _check_arguments():
    """Exactly one input mode, with the matching output argument."""
    if (par["input"] is None) == (par["input_table"] is None):
        raise ValueError(
            "Exactly one of --input (h5mu) and --input_table (CSV) must be given."
        )
    if par["input"] is not None:
        if par["output"] is None:
            raise ValueError("--output is required when --input is given.")
        if par["metadata"] is not None:
            raise ValueError(
                "--metadata only applies to --input_table; with --input the labels "
                "are read from .obs."
            )
    else:
        if par["output_table"] is None:
            raise ValueError(
                "--output_table is required when --input_table is given."
            )


def _label_source(par):
    """Where the label columns come from, for error messages."""
    return ".obs" if par["input"] is not None else "--metadata"


def _resolve_start_group(adata, par):
    """Return a single observation identifier to use as the trajectory root."""
    if par["start_group"] is not None and par["start_group_cluster"] is not None:
        raise ValueError(
            "--start_group and --start_group_cluster are mutually exclusive; "
            f"got --start_group '{par['start_group']}' and --start_group_cluster "
            f"'{par['start_group_cluster']}'. Provide exactly one."
        )
    if par["start_group_cluster"] is not None:
        obs_key = par["start_group_column"]
        cluster = par["start_group_cluster"]
        if obs_key not in adata.obs.columns:
            raise ValueError(
                f"--start_group_column '{obs_key}' not found in {_label_source(par)}. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        if not (adata.obs[obs_key] == cluster).any():
            raise ValueError(
                f"No observations found for '{cluster}' in "
                f"{_label_source(par)} column '{obs_key}'."
            )
        # early_cell(ad, celltype, celltype_column) - requires DM_EigenVectors_multiscaled
        start_group = palantir.utils.early_cell(adata, cluster, celltype_column=obs_key)
        logger.info(
            "Auto-selected root '%s' from population '%s' (column '%s')",
            start_group,
            cluster,
            obs_key,
        )
        return start_group
    elif par["start_group"] is not None:
        if par["start_group"] not in adata.obs_names:
            raise ValueError(
                f"--start_group '{par['start_group']}' is not one of the "
                f"{adata.n_obs} observation identifiers."
            )
        return par["start_group"]
    else:
        raise ValueError(
            "Either --start_group or --start_group_cluster must be provided."
        )


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
                f"{_label_source(par)}. Available columns: {list(adata.obs.columns)}"
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


def _read_h5mu_input():
    """Return (adata, obsm_key) for h5mu mode."""
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


def _read_table_input():
    """Build an in-memory AnnData from the embedding table (+ optional metadata)."""
    logger.info("Reading table %s", par["input_table"])
    embedding = read_group_table(par["input_table"], par["id_column"], "--input_table")
    id_column = embedding.index.name

    obs = pd.DataFrame(index=embedding.index)
    if par["metadata"] is not None:
        logger.info("Reading metadata %s", par["metadata"])
        meta_df = pd.read_csv(par["metadata"])
        if id_column not in meta_df.columns:
            raise ValueError(
                f"Identifier column '{id_column}' not found in --metadata "
                f"'{par['metadata']}'. Available: {list(meta_df.columns)}"
            )
        meta_df[id_column] = meta_df[id_column].astype(str)
        meta_df = meta_df.drop_duplicates(subset=[id_column]).set_index(id_column)
        missing = embedding.index.difference(meta_df.index)
        if len(missing) > 0:
            raise ValueError(
                f"--metadata '{par['metadata']}' is missing {len(missing)} "
                f"identifier(s) present in --input_table, e.g. {list(missing[:5])}."
            )
        obs = meta_df.loc[embedding.index]

    adata = ad.AnnData(
        X=np.zeros((embedding.shape[0], 0), dtype="float32"),
        obs=obs,
        obsm={TABLE_EMBEDDING_KEY: embedding.to_numpy(dtype=float)},
    )
    adata.uns["_id_column"] = id_column
    return adata, TABLE_EMBEDDING_KEY


def _write_table_output(adata):
    """Flatten the Palantir results back into one CSV."""
    out = pd.DataFrame(index=adata.obs_names)
    out[par["pseudotime_column"]] = adata.obs[par["pseudotime_column"]].to_numpy()
    out[par["entropy_column"]] = adata.obs[par["entropy_column"]].to_numpy()

    fate = adata.obsm.get(par["obsm_fate_probabilities"])
    if fate is not None:
        fate_df = pd.DataFrame(fate)
        fate_df.index = adata.obs_names
        for col in fate_df.columns:
            out[f"fate_{col}"] = fate_df[col].to_numpy()

    waypoints = adata.uns.get(par["uns_waypoints"])
    if waypoints is not None:
        out["palantir_waypoint"] = adata.obs_names.isin(np.asarray(waypoints))

    out.index = pd.Index(out.index, name=adata.uns["_id_column"])
    write_group_table(out, par["output_table"])
    logger.info(
        "Written Palantir results to %s (shape %s)", par["output_table"], out.shape
    )


def main():
    _check_arguments()
    np.random.seed(par["seed"])

    if par["input"] is not None:
        adata, obsm_key = _read_h5mu_input()
    else:
        adata, obsm_key = _read_table_input()

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
    start_group = _resolve_start_group(adata, par)
    logger.info("Trajectory root: %s", start_group)

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
        early_cell=start_group,
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
    if par["input"] is not None:
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
    else:
        _write_table_output(adata)


if __name__ == "__main__":
    main()
