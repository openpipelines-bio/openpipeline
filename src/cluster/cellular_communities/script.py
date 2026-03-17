import sys
import numpy as np
import pandas as pd
import mudata as mu
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

## VIASH START
par = {
    "input": "test_data/communities_input.h5mu",
    "modality": "rna",
    "obs_subpopulation": "subpopulation",
    "uns_proportions": "proportions",
    "uns_dynamics": "dynamics",
    "n_communities": 3,
    "alpha": 0.5,
    "method": "hierarchical",
    "linkage": "ward",
    "output": "test_data/communities_output.h5mu",
    "obs_community_id": "community_id",
    "uns_output": "cellular_communities",
    "output_compression": None,
}
meta = {
    "resources_dir": "src/cluster/cellular_communities/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _cooccurrence_similarity(prop_df):
    """Pearson correlation matrix of participant proportion vectors (subpop × subpop)."""
    # prop_df: participants × subpopulations
    corr = prop_df.corr(method="pearson")
    # clip to [-1, 1] and fill NaN (zero-variance subpops) with 0
    corr = corr.clip(-1, 1).fillna(0)
    return corr


def _dynamics_similarity(dynamics, subpops):
    """Pearson correlation of fitted proportion curves (subpop × subpop)."""
    n = len(subpops)
    sim = np.eye(n)
    curves = {}
    for sp in subpops:
        if sp in dynamics:
            curves[sp] = np.array(dynamics[sp]["proportion_fitted"])

    for i, sp_i in enumerate(subpops):
        for j, sp_j in enumerate(subpops):
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
    return pd.DataFrame(sim, index=subpops, columns=subpops)


def _cluster_hierarchical(dist_mat, n_communities, link_method):
    """Ward linkage hierarchical clustering; returns integer cluster labels."""
    condensed = squareform(dist_mat, checks=False)
    condensed = np.clip(condensed, 0, None)   # numerical noise may give tiny negatives
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
    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    # ── validate ─────────────────────────────────────────────────────────────
    subpop_col = par["obs_subpopulation"]
    if subpop_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{subpop_col}' not found in .obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    for key in (par["uns_proportions"], par["uns_dynamics"]):
        if key not in adata.uns:
            raise ValueError(
                f"Key '{key}' not found in .uns. "
                f"Available: {list(adata.uns.keys())}"
            )

    # ── load data ─────────────────────────────────────────────────────────────
    prop_df = pd.DataFrame(adata.uns[par["uns_proportions"]])  # participants × subpops
    dynamics = adata.uns[par["uns_dynamics"]]
    subpops = prop_df.columns.tolist()
    logger.info(
        "Proportion matrix: %d participants × %d subpopulations.",
        *prop_df.shape,
    )

    # ── similarity matrices ───────────────────────────────────────────────────
    logger.info("Computing co-occurrence similarity (Pearson correlation).")
    co_sim = _cooccurrence_similarity(prop_df)

    logger.info("Computing dynamics similarity (curve correlation).")
    dyn_sim = _dynamics_similarity(dynamics, subpops)

    alpha = par["alpha"]
    combined = alpha * co_sim.values + (1.0 - alpha) * dyn_sim.values
    combined = np.clip(combined, -1, 1)
    combined_df = pd.DataFrame(combined, index=subpops, columns=subpops)
    logger.info("Combined similarity matrix (alpha=%.2f) computed.", alpha)

    # ── clustering ───────────────────────────────────────────────────────────
    dist_mat = 1.0 - combined_df
    n_comm = par["n_communities"]
    method = par["method"]

    if method == "hierarchical":
        logger.info(
            "Hierarchical clustering (linkage=%s, n_communities=%d).",
            par["linkage"], n_comm,
        )
        community_labels = _cluster_hierarchical(dist_mat.values, n_comm, par["linkage"])
    else:
        logger.info("Spectral clustering (n_communities=%d).", n_comm)
        community_labels = _cluster_spectral(combined_df, n_comm)

    subpop_to_community = dict(zip(subpops, community_labels))
    logger.info("Community assignments: %s", subpop_to_community)

    # ── assign to cells ───────────────────────────────────────────────────────
    comm_col = par["obs_community_id"]
    adata.obs[comm_col] = adata.obs[subpop_col].map(subpop_to_community).fillna("unassigned")
    logger.info("Assigned community labels to .obs['%s'].", comm_col)

    # ── store metadata in uns ─────────────────────────────────────────────────
    uns_key = par["uns_output"]
    adata.uns[uns_key] = {
        "subpopulation_communities": subpop_to_community,
        "n_communities": n_comm,
        "alpha": alpha,
        "method": method,
        "co_occurrence_similarity": co_sim.to_dict(),
        "dynamics_similarity": dyn_sim.to_dict(),
        "combined_similarity": combined_df.to_dict(),
    }
    logger.info("Stored community metadata in .uns['%s'].", uns_key)

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
