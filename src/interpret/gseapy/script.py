import sys
import os
import pandas as pd
import gseapy as gp

## VIASH START
par = {
    "input": "deseq2_results.csv",
    "gene_column": None,
    "fc_column": "log2FoldChange",
    "pval_column": "padj",
    "gene_sets_file": None,
    "method": "prerank",
    "gene_sets": ["MSigDB_Hallmark_2020"],
    "pval_threshold": 0.05,
    "fc_threshold": 0.0,
    "organism": "Human",
    "min_size": 15,
    "max_size": 500,
    "permutation_num": 1000,
    "seed": 42,
    "output": "pathway_enrichment.csv",
}
meta = {
    "cpus": 2,
    "resources_dir": "src/interpret/gseapy/",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def _collect_gene_sets(par):
    """Return the list of gene set libraries: Enrichr names plus staged GMT files.

    GMT files arrive through `--gene_sets_file` (`type: file`), so viash stages them for
    both the executable and the nextflow runner. `--gene_sets` only ever holds Enrichr
    library names, which gseapy resolves over the network.
    """
    gene_sets = list(par["gene_sets"] or []) + list(par["gene_sets_file"] or [])
    if not gene_sets:
        raise ValueError(
            "No gene sets provided. Give at least one Enrichr library name via "
            "--gene_sets or at least one GMT file via --gene_sets_file."
        )
    return gene_sets


def _check_numeric_column(de, column, argument):
    """Fail early when a required DE column is missing or not numeric."""
    if column not in de.columns:
        raise ValueError(
            f"{argument} '{column}' not found in --input. "
            f"Available columns: {list(de.columns)}"
        )
    if not pd.api.types.is_numeric_dtype(de[column]):
        raise ValueError(
            f"{argument} '{column}' is not numeric (dtype: {de[column].dtype}). "
            "Check the DE table for header or separator problems."
        )
    if de[column].isna().all():
        raise ValueError(f"{argument} '{column}' contains only NA values.")


def _load_de_table(csv_path, gene_column):
    """Load DESeq2 CSV and return a DataFrame indexed by gene name."""
    de = pd.read_csv(csv_path, index_col=0)
    if de.empty:
        raise ValueError(f"--input '{csv_path}' contains no rows.")
    if gene_column:
        if gene_column not in de.columns:
            raise ValueError(
                f"--gene_column '{gene_column}' not found in DE CSV. "
                f"Available columns: {list(de.columns)}"
            )
        de = de.set_index(gene_column)
    n_duplicated = de.index.duplicated().sum()
    if n_duplicated:
        raise ValueError(
            f"{n_duplicated} duplicated gene name(s) in --input, "
            "for example: "
            f"{list(de.index[de.index.duplicated()][:5])}"
        )
    return de


def _run_prerank(de, gene_sets, par, n_jobs):
    """Run pre-ranked GSEA for each gene set library."""
    results = {}
    ranking = de[par["fc_column"]].dropna().sort_values(ascending=False)
    if ranking.empty:
        raise ValueError(f"Ranking column '{par['fc_column']}' has no non-NA values.")
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
        df = res.res2d.copy()
        logger.info("  %d terms for %s", len(df), label)
        results[label] = df
    return results


def _run_ora(de, gene_sets, par, n_jobs):
    """Run over-representation analysis for each gene set library."""
    results = {}
    sig_mask = de[par["pval_column"]] < par["pval_threshold"]
    if par["fc_threshold"] > 0:
        sig_mask = sig_mask & (de[par["fc_column"]].abs() >= par["fc_threshold"])
    gene_list = de.index[sig_mask & de[par["pval_column"]].notna()].tolist()
    if not gene_list:
        raise ValueError(
            f"No significant genes with --pval_column '{par['pval_column']}' < "
            f"{par['pval_threshold']} and |{par['fc_column']}| >= {par['fc_threshold']}. "
            "ORA has nothing to test; loosen the thresholds or check the DE table."
        )
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
        df = res.res2d.copy()
        logger.info("  %d terms for %s", len(df), label)
        results[label] = df
    return results


def main():
    logger.info("Reading DE results from %s", par["input"])
    de = _load_de_table(par["input"], par["gene_column"])
    logger.info("  %d genes loaded", len(de))

    _check_numeric_column(de, par["fc_column"], "--fc_column")

    gene_sets = _collect_gene_sets(par)
    n_jobs = max(1, (meta.get("cpus") or 1))

    if par["method"] == "prerank":
        enrichment_results = _run_prerank(de, gene_sets, par, n_jobs)
    elif par["method"] == "ora":
        _check_numeric_column(de, par["pval_column"], "--pval_column")
        enrichment_results = _run_ora(de, gene_sets, par, n_jobs)
    else:
        raise ValueError(
            f"Unknown method '{par['method']}'. Choose 'prerank' or 'ora'."
        )

    # One long table for all libraries; the library is a column, not a file name
    combined = pd.concat(
        [
            df.assign(gene_set_library=label, method=par["method"])
            for label, df in enrichment_results.items()
        ],
        ignore_index=True,
    )
    lead_cols = ["gene_set_library", "method"]
    combined = combined[lead_cols + [c for c in combined.columns if c not in lead_cols]]

    combined.to_csv(par["output"], index=False)
    logger.info(
        "Written %d rows for %d gene set librar(y/ies) to %s",
        len(combined),
        len(enrichment_results),
        par["output"],
    )


if __name__ == "__main__":
    main()
