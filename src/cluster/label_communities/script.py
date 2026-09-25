import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

## VIASH START
par = {
    "input": "proportions.csv",
    "dynamics": "dynamics.csv",
    "id_column": None,
    "n_communities": 3,
    "alpha": 0.5,
    "correlation_method": "spearman",
    "method": "hierarchical",
    "linkage": "ward",
    "output": "communities.csv",
    "output_similarity": None,
}
meta = {
    "resources_dir": "src/cluster/label_communities/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from group_table import read_group_table

logger = setup_logger()


def _cooccurrence_similarity(prop_df, method):
    """Correlation matrix of group proportion vectors (label x label)."""
    # prop_df: groups x labels; a zero-variance label gives NaN
    corr = prop_df.corr(method=method)
    corr = corr.fillna(0)
    return corr


def _read_dynamics(path, labels):
    """Long-format curve table -> {label: fitted curve}, ordered by pseudotime."""
    df = pd.read_csv(path)
    required = {"label", "pseudotime", "proportion_fitted"}
    missing = required.difference(df.columns)
    df = pd.read_csv(path)
    required = {"label", "pseudotime", "proportion_fitted"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"--dynamics '{path}' is missing column(s) {sorted(missing)}. "
            f"Available: {list(df.columns)}"
        )
    
    curves = {}
    for label, group in df.groupby("label"):
        curves[str(label)] = group.sort_values("pseudotime")[
            "proportion_fitted"
        ].to_numpy(dtype=float)
    
    absent = [lab for lab in labels if lab not in curves]
    if absent:
        logger.warning(
            "%d label(s) have no fitted curve in --dynamics and get zero dynamics "
            "similarity: %s",
            len(absent),
            absent[:10],
        )
    
    return curves


def _dynamics_similarity(curves, labels):
    """Pearson correlation of fitted proportion curves (label x label)."""
    n = len(labels)
    sim = np.eye(n)
    for i, sp_i in enumerate(labels):
        for j, sp_j in enumerate(labels):
            if i >= j:
                continue
            if sp_i in curves and sp_j in curves:
                ci, cj = curves[sp_i], curves[sp_j]
                min_len = min(len(ci), len(cj))
                if min_len < 3:
                    r = 0.0
                else:
                    r, _ = pearsonr(ci[:min_len], cj[:min_len])
                    r = float(np.clip(r, -1, 1)) if np.isfinite(r) else 0.0
            else:
                r = 0.0
            sim[i, j] = r
            sim[j, i] = r
    return pd.DataFrame(sim, index=labels, columns=labels)


def _cluster_hierarchical(dist_mat, n_communities, link_method):
    """Hierarchical clustering; returns integer cluster labels."""
    condensed = squareform(dist_mat, checks=False)
    condensed = np.clip(condensed, 0, None)  # numerical noise may give tiny negatives
    Z = linkage(condensed, method=link_method)
    labels = fcluster(Z, t=n_communities, criterion="maxclust")
    return labels.astype(str)


def _cluster_spectral(sim_mat, n_communities):
    """Spectral clustering on a similarity matrix."""
    from sklearn.cluster import SpectralClustering

    sc = SpectralClustering(
        n_clusters=n_communities,
        affinity="precomputed",
        random_state=42,
        n_init=10,
    )
    arr = np.clip(sim_mat.values, 0, None)
    labels = sc.fit_predict(arr)
    return (labels + 1).astype(str)


def _write_similarity(path, co_sim, dyn_sim, combined_df, labels):
    rows = []
    for i, sp_i in enumerate(labels):
        for sp_j in labels[i + 1 :]:
            rows.append(
                {
                    "label_1": sp_i,
                    "label_2": sp_j,
                    "co_occurrence": float(co_sim.loc[sp_i, sp_j]),
                    "dynamics": float(dyn_sim.loc[sp_i, sp_j]),
                    "combined": float(combined_df.loc[sp_i, sp_j]),
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)
    logger.info("Written %d label pairs to %s.", len(rows), path)


def main():
    alpha = par["alpha"]
    if alpha < 1.0 and par["dynamics"] is None:
        raise ValueError(
            f"--dynamics is required unless --alpha is 1.0; got --alpha {alpha}."
        )

    logger.info("Reading proportions from %s", par["input"])
    prop_df = read_group_table(par["input"], par["id_column"], "--input")
    labels = [str(c) for c in prop_df.columns]
    prop_df.columns = labels
    logger.info("Proportion matrix: %d groups x %d labels.", *prop_df.shape)

    if par["n_communities"] > len(labels):
        raise ValueError(
            f"--n_communities ({par['n_communities']}) exceeds the number of labels "
            f"in --input ({len(labels)})."
        )

    # -- similarity matrices ---------------------------------------------------
    logger.info(
        "Computing co-occurrence similarity (%s correlation).",
        par["correlation_method"],
    )
    co_sim = _cooccurrence_similarity(prop_df, par["correlation_method"])

    if par["dynamics"] is not None:
        logger.info("Computing dynamics similarity from %s.", par["dynamics"])
        curves = _read_dynamics(par["dynamics"], labels)
        dyn_sim = _dynamics_similarity(curves, labels)
    else:
        logger.info("No --dynamics given; dynamics similarity is the identity.")
        dyn_sim = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)

    combined = alpha * co_sim.values + (1.0 - alpha) * dyn_sim.values
    combined = np.clip(combined, -1, 1)
    combined_df = pd.DataFrame(combined, index=labels, columns=labels)
    logger.info("Combined similarity matrix (alpha=%.2f) computed.", alpha)

    # -- clustering -----------------------------------------------------------
    dist_mat = 1.0 - combined_df
    n_comm = par["n_communities"]
    method = par["method"]

    if method == "hierarchical":
        logger.info(
            "Hierarchical clustering (linkage=%s, n_communities=%d).",
            par["linkage"],
            n_comm,
        )
        community_labels = _cluster_hierarchical(
            dist_mat.values, n_comm, par["linkage"]
        )
    else:
        logger.info("Spectral clustering (n_communities=%d).", n_comm)
        community_labels = _cluster_spectral(combined_df, n_comm)

    out = pd.DataFrame({"label": labels, "community_id": community_labels})
    out.to_csv(par["output"], index=False)
    logger.info(
        "Assigned %d labels to %d communities; written to %s.",
        len(labels),
        out["community_id"].nunique(),
        par["output"],
    )

    if par["output_similarity"]:
        _write_similarity(
            par["output_similarity"], co_sim, dyn_sim, combined_df, labels
        )


if __name__ == "__main__":
    main()
