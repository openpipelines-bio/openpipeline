import sys
import numpy as np
import mudata as mu

## VIASH START
par = {
    "input": "test_data/via_input.h5mu",
    "modality": "rna",
    "obsm_key": "X_phate",
    "obs_cluster": "leiden",
    "root_user": "0",
    "knn": 15,
    "random_seed": 42,
    "n_iter_leiden": 5,
    "output": "test_data/via_output.h5mu",
    "obs_pseudotime": "via_pseudotime",
    "uns_graph": "via_graph",
    "output_compression": None,
}
meta = {
    "resources_dir": "src/trajectory/via/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _parse_root(root_str, cluster_labels):
    """Convert root_user string to the format VIA expects.

    Returns a list with either a cluster label string or an integer cell index.
    """
    try:
        return [int(root_str)]
    except ValueError:
        if root_str not in cluster_labels:
            raise ValueError(
                f"root_user '{root_str}' is not a valid cell index (integer) "
                f"nor a cluster label. Available clusters: {sorted(set(cluster_labels))}"
            )
        return [root_str]


def main():
    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    # ── validate ──────────────────────────────────────────────────────────────
    obsm_key = par["obsm_key"]
    if obsm_key not in adata.obsm:
        raise ValueError(
            f"obsm key '{obsm_key}' not found. Available: {list(adata.obsm.keys())}"
        )
    obs_cluster = par["obs_cluster"]
    if obs_cluster not in adata.obs.columns:
        raise ValueError(
            f"obs column '{obs_cluster}' not found. Available: {list(adata.obs.columns)}"
        )

    # ── prepare inputs ────────────────────────────────────────────────────────
    embedding = np.array(adata.obsm[obsm_key], dtype=float)
    cluster_labels = adata.obs[obs_cluster].astype(str).tolist()
    root_user = _parse_root(par["root_user"], cluster_labels)

    logger.info(
        "Running VIA on %d cells, embedding %s, root=%s.",
        len(cluster_labels), embedding.shape, root_user,
    )

    # ── run VIA ───────────────────────────────────────────────────────────────
    import pyVIA.core as via_core
    import pyVIA.utils_via as via_utils
    from scipy.sparse import csr_matrix as _csr

    def _patched_get_sparse_from_igraph(graph, weight_attr=None):
        """pyVIA 0.2.4 bug-fix: zip(*[]) fails with newer scipy when edge list is empty."""
        n = graph.vcount()
        edges = list(graph.get_edgelist())
        if not edges:
            return _csr((n, n))
        weights = graph.es[weight_attr] if weight_attr else [1.0] * len(edges)
        rows, cols = zip(*edges)
        return _csr((weights, (list(rows), list(cols))), shape=(n, n))

    # Patch in both module namespaces: core (from utils_via import *) and utils_via itself
    via_core.get_sparse_from_igraph = _patched_get_sparse_from_igraph
    via_utils.get_sparse_from_igraph = _patched_get_sparse_from_igraph

    v = via_core.VIA(
        data=embedding,
        true_label=cluster_labels,
        root_user=root_user,
        knn=par["knn"],
        random_seed=par["random_seed"],
        n_iter_leiden=par["n_iter_leiden"],
        num_threads=meta.get("cpus") or 1,
        jac_weighted_edges=False,   # disable Jaccard pruning that can empty the graph
        keep_all_local_dist=True,   # retain all local edges
    )
    v.run_VIA()
    logger.info("VIA run complete.")

    # ── extract pseudotime ────────────────────────────────────────────────────
    pt = np.array(v.single_cell_pt_markov, dtype=float)
    adata.obs[par["obs_pseudotime"]] = pt
    logger.info(
        "Pseudotime range: [%.4f, %.4f].", float(pt.min()), float(pt.max())
    )

    # ── extract graph ─────────────────────────────────────────────────────────
    # v.edgelist: list of (source_cluster, target_cluster, weight) tuples
    edge_list = getattr(v, "edgelist", None) or []
    adata.uns[par["uns_graph"]] = {
        "source": [int(e[0]) for e in edge_list],
        "target": [int(e[1]) for e in edge_list],
        "weight": [float(e[2]) for e in edge_list],
    }
    logger.info("Stored VIA graph with %d edges in .uns['%s'].", len(edge_list), par["uns_graph"])

    # ── write output ──────────────────────────────────────────────────────────
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
