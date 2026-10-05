import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

## VIASH START
par = {
    "input": ["shift_batch_0001.csv", "shift_batch_0002.csv"],
    "shift_column": "shift_healthy",
    "min_coverage": 350,
    "n_random": 1000,
    "seed": 41,
    "pvalue_cutoff": 0.05,
    "top_n": None,
    "output_ranked": "ranked_genes.csv",
    "output_top": "top_genes.csv",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()

OUTPUT_COLUMNS = ["rank", "gene_id", "gene_name", "median_shift", "pvalue", "n_cells"]


def write(ranked, n_top):
    ranked.to_csv(par["output_ranked"], index=False)
    ranked.iloc[:n_top].to_csv(par["output_top"], index=False)


def main():
    shifts = pd.concat([pd.read_csv(path) for path in par["input"]], ignore_index=True)
    logger.info("Read %i files, %i knockouts", len(par["input"]), shifts.shape[0])

    required = ["cell_id", "gene_id", "gene_name", par["shift_column"]]
    missing = [column for column in required if column not in shifts.columns]
    if missing:
        raise ValueError(
            f"Columns {missing} not found in the similarity shift input. "
            f"Available: {list(shifts.columns)}"
        )
    shifts["gene_id"] = shifts["gene_id"].astype(str)

    by_gene = shifts.groupby("gene_id", sort=False)
    genes = pd.DataFrame(
        {
            "gene_name": by_gene["gene_name"].first(),
            "median_shift": by_gene[par["shift_column"]].median(),
            "n_cells": by_gene["cell_id"].nunique(),
        }
    )
    covered = genes[genes["n_cells"] > par["min_coverage"]]
    logger.info(
        "%i of %i genes are knocked out in more than %i cells",
        covered.shape[0],
        genes.shape[0],
        par["min_coverage"],
    )
    if covered.empty:
        logger.warning("No gene passes the coverage filter, writing empty tables")
        write(pd.DataFrame(columns=OUTPUT_COLUMNS), 0)
        return

    rng = np.random.default_rng(par["seed"])
    n_random = min(par["n_random"], shifts.shape[0])
    baseline = shifts[par["shift_column"]].to_numpy()[
        rng.permutation(shifts.shape[0])[:n_random]
    ]
    logger.info("Random baseline of %i shifts (seed %i)", n_random, par["seed"])

    values = {gene: group[par["shift_column"]].to_numpy() for gene, group in by_gene}
    ranked = (
        covered.sort_values("median_shift", ascending=False, kind="stable")
        .rename_axis("gene_id")
        .reset_index()
    )
    ranked["pvalue"] = [
        mannwhitneyu(values[gene], baseline, alternative="greater").pvalue
        for gene in ranked["gene_id"]
    ]
    ranked["rank"] = np.arange(ranked.shape[0]) + 1
    ranked = ranked[OUTPUT_COLUMNS]

    if par["top_n"] is not None:
        n_top = min(par["top_n"], ranked.shape[0])
    else:
        significant = (ranked["pvalue"] < par["pvalue_cutoff"]).to_numpy()
        # length of the uninterrupted run of significant genes from rank 1
        n_top = ranked.shape[0] if significant.all() else int(np.argmin(significant))
    logger.info("Ranked %i genes, %i in the top table", ranked.shape[0], n_top)
    write(ranked, n_top)


if __name__ == "__main__":
    main()
