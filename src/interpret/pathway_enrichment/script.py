import sys
import os
import pandas as pd
import mudata as mu
import gseapy as gp

## VIASH START
par = {
    "input": "src/interpret/pathway_enrichment/test_data/input.h5mu",
    "modality": "rna",
    "input_degenes": "src/interpret/pathway_enrichment/test_data/deseq2_results.csv",
    "gene_column": None,
    "method": "prerank",
    "gene_sets": ["MSigDB_Hallmark_2020"],
    "fc_column": "log2FoldChange",
    "pval_column": "padj",
    "pval_threshold": 0.05,
    "fc_threshold": 0.0,
    "organism": "Human",
    "min_size": 15,
    "max_size": 500,
    "permutation_num": 1000,
    "seed": 42,
    "uns_key": "pathway_enrichment",
    "output": "output.h5mu",
    "output_csv_dir": "pathway_enrichment_results/",
    "output_compression": None,
}
meta = {
    "cpus": 2,
    "resources_dir": "src/interpret/pathway_enrichment/",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def _resolve_gene_set(gs, input_paths=None):
    """Return the accessible path for a GMT file, or gs unchanged for Enrichr names.

    Inside a Docker container viash mounts the host filesystem under /viash_automount.
    File-type arguments (--input, --input_degenes) are automatically converted to their
    /viash_automount/... paths, but string-type arguments like --gene_sets are not.

    Strategy:
    1. If gs is already an accessible file path, return it as-is.
    2. If absolute: prepend /viash_automount.
    3. If relative: walk up the directory tree of each known automounted input path
       until the relative path resolves to an existing file (finds the host CWD).
    4. Fall through: treat as an Enrichr library name.
    """
    if os.path.isfile(gs):
        return gs

    mount_root = "/viash_automount"
    if not os.path.isdir(mount_root):
        return gs  # not running inside Docker

    if os.path.isabs(gs):
        candidate = mount_root + gs
        if os.path.isfile(candidate):
            return candidate
        return gs

    # Relative path: infer host CWD by walking up from automounted input paths
    for inp in (input_paths or []):
        if not (inp and inp.startswith(mount_root)):
            continue
        d = os.path.dirname(inp)
        while d and d != mount_root:
            candidate = os.path.join(d, gs)
            if os.path.isfile(candidate):
                return candidate
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent

    return gs  # treat as Enrichr library name


def _load_de_table(csv_path, gene_column):
    """Load DESeq2 CSV and return a DataFrame indexed by gene name."""
    de = pd.read_csv(csv_path, index_col=0)
    if gene_column:
        if gene_column not in de.columns:
            raise ValueError(
                f"--gene_column '{gene_column}' not found in DE CSV. "
                f"Available columns: {list(de.columns)}"
            )
        de = de.set_index(gene_column)
    return de


def _run_prerank(de, gene_sets, output_dir, par, n_jobs):
    """Run pre-ranked GSEA for each gene set library."""
    results = {}
    ranking = (
        de[par["fc_column"]]
        .dropna()
        .sort_values(ascending=False)
    )
    if ranking.empty:
        raise ValueError(
            f"Ranking column '{par['fc_column']}' has no non-NA values."
        )
    for gs in gene_sets:
        logger.info("prerank GSEA with gene set: %s", gs)
        label = os.path.splitext(os.path.basename(gs))[0] if os.path.isfile(gs) else gs
        res = gp.prerank(
            rnk=ranking,
            gene_sets=gs,
            min_size=par["min_size"],
            max_size=par["max_size"],
            permutation_num=par["permutation_num"],
            seed=par["seed"],
            threads=n_jobs,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        df = res.res2d
        out_csv = os.path.join(output_dir, f"{label}.prerank.csv")
        df.to_csv(out_csv)
        logger.info("  saved %d terms to %s", len(df), out_csv)
        results[label] = df.to_dict(orient="list")
    return results


def _run_ora(de, gene_sets, output_dir, par, n_jobs):
    """Run over-representation analysis for each gene set library."""
    results = {}
    sig_mask = de[par["pval_column"]] < par["pval_threshold"]
    if par["fc_threshold"] > 0:
        sig_mask = sig_mask & (de[par["fc_column"]].abs() >= par["fc_threshold"])
    gene_list = de.index[sig_mask & de[par["pval_column"]].notna()].tolist()
    if not gene_list:
        logger.warning(
            "No significant genes found with pval_threshold=%.3f, fc_threshold=%.2f. "
            "Skipping ORA.",
            par["pval_threshold"],
            par["fc_threshold"],
        )
        return results
    logger.info("ORA with %d significant genes", len(gene_list))
    for gs in gene_sets:
        label = os.path.splitext(os.path.basename(gs))[0] if os.path.isfile(gs) else gs
        logger.info("ORA with gene set: %s", gs)
        res = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gs,
            organism=par["organism"].lower(),
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        df = res.res2d
        out_csv = os.path.join(output_dir, f"{label}.ora.csv")
        df.to_csv(out_csv)
        logger.info("  saved %d terms to %s", len(df), out_csv)
        results[label] = df.to_dict(orient="list")
    return results


def main():
    logger.info("Reading DE results from %s", par["input_degenes"])
    de = _load_de_table(par["input_degenes"], par["gene_column"])
    logger.info("  %d genes loaded", len(de))

    if par["fc_column"] not in de.columns:
        raise ValueError(
            f"--fc_column '{par['fc_column']}' not found. "
            f"Available: {list(de.columns)}"
        )

    os.makedirs(par["output_csv_dir"], exist_ok=True)

    n_jobs = max(1, (meta.get("cpus") or 1))
    input_paths = [par.get("input"), par.get("input_degenes")]
    gene_sets = [_resolve_gene_set(gs, input_paths=input_paths) for gs in par["gene_sets"]]

    if par["method"] == "prerank":
        enrichment_results = _run_prerank(de, gene_sets, par["output_csv_dir"], par, n_jobs)
    elif par["method"] == "ora":
        if par["pval_column"] not in de.columns:
            raise ValueError(
                f"--pval_column '{par['pval_column']}' not found. "
                f"Available: {list(de.columns)}"
            )
        enrichment_results = _run_ora(de, gene_sets, par["output_csv_dir"], par, n_jobs)
    else:
        raise ValueError(f"Unknown method '{par['method']}'. Choose 'prerank' or 'ora'.")

    logger.info("Reading input h5mu from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    mod = mdata.mod[par["modality"]]

    if par["uns_key"] not in mod.uns:
        mod.uns[par["uns_key"]] = {}
    mod.uns[par["uns_key"]].update(enrichment_results)
    logger.info(
        "Stored enrichment results for %d gene set(s) in .uns['%s']",
        len(enrichment_results),
        par["uns_key"],
    )

    logger.info("Writing output to %s", par["output"])
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
