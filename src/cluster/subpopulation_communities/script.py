import sys
import numpy as np
import pandas as pd
import mudata as mu
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

## VIASH START
par = {
    "input": "proportions.h5mu",
    "modality": "proportions",
    "varm_dynamics": "proportion_dynamics",
    "n_communities": 3,
    "alpha": 0.5,
    "correlation_method": "spearman",
    "method": "hierarchical",
    "linkage": "ward",
    "output": "communities.h5mu",
    "var_output": "community_id",
    "varp_prefix": "community_",
    "output_compression": None,
}
meta = {
    "resources_dir": "src/cluster/subpopulation_communities/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _cooccurrence_similarity(prop_df, method):
    """Correlation matrix of group proportion vectors (label x label)."""
    # prop_df: groups x labels; a zero-variance label gives NaN
    corr = prop_df.corr(method=method)
    corr = corr.fillna(0)
    return corr


def _read_dynamics(adata, key):
    """Fitted curves in .varm -> {label: fitted curve}; labels without a fit are left out."""
    if key not in adata.varm:
        raise ValueError(
            f"--varm_dynamics '{key}' not found in .varm, but it is required unless "
            f"--alpha is 1.0. Available keys: {list(adata.varm.keys())}"
        )
    values = np.asarray(adata.varm[key], dtype=float)

    curves = {}
    absent = []
    for label, curve in zip(adata.var_names, values):
        if np.isfinite(curve).all():
            curves[str(label)] = curve
        else:
            absent.append(str(label))
    if absent:
        logger.warning(
            "%d label(s) have no fitted curve in .varm['%s'] and get zero dynamics "
            "similarity: %s",
            len(absent),
            key,
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


def main():
    alpha = par["alpha"]

    logger.info("Reading %s, modality '%s'.", par["input"], par["modality"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])
    labels = [str(c) for c in adata.var_names]
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    prop_df = pd.DataFrame(X, index=adata.obs_names, columns=labels, dtype=float)
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

    if alpha < 1.0:
        logger.info(
            "Computing dynamics similarity from .varm['%s'].", par["varm_dynamics"]
        )
        curves = _read_dynamics(adata, par["varm_dynamics"])
        dyn_sim = _dynamics_similarity(curves, labels)
    else:
        logger.info("--alpha is 1.0; dynamics similarity is not used.")
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

    adata.var[par["var_output"]] = pd.Categorical(community_labels)
    prefix = par["varp_prefix"]
    adata.varp[f"{prefix}co_occurrence"] = co_sim.to_numpy(dtype=float)
    adata.varp[f"{prefix}dynamics"] = dyn_sim.to_numpy(dtype=float)
    adata.varp[f"{prefix}combined"] = combined_df.to_numpy(dtype=float)
    logger.info(
        "Assigned %d labels to %d communities in .var['%s'].",
        len(labels),
        len(set(community_labels)),
        par["var_output"],
    )

    logger.info("Writing output to %s.", par["output"])
    write_h5ad_to_h5mu_with_compression(
        par["output"],
        par["input"],
        par["modality"],
        adata,
        par["output_compression"],
    )


if __name__ == "__main__":
    main()
