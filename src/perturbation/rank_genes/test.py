import re
import sys
from subprocess import CalledProcessError

import numpy as np
import pandas as pd
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/perturbation/rank_genes/rank_genes",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/rank_genes/config.vsh.yaml",
}
## VIASH END

N_CELLS = 40


@pytest.fixture
def shift_paths(tmp_path):
    """Two batches of 20 cells. Gene UP moves every cell towards healthy,
    ZERO does nothing, DOWN moves away, and RARE is only knocked out in 3 cells."""
    rng = np.random.default_rng(3)
    rows = []
    for cell in range(N_CELLS):
        for gene, centre in [("UP", 0.05), ("ZERO", 0.0), ("DOWN", -0.05)]:
            rows.append(
                (f"cell_{cell}", f"ENSG_{gene}", gene, centre + rng.normal(0, 0.005))
            )
        if cell < 3:
            rows.append((f"cell_{cell}", "ENSG_RARE", "RARE", 0.5))
    shifts = pd.DataFrame(
        rows, columns=["cell_id", "gene_id", "gene_name", "shift_healthy"]
    )
    shifts["shift_disease"] = -shifts["shift_healthy"]
    paths = []
    for batch, cells in enumerate([range(0, 20), range(20, 40)]):
        path = tmp_path / f"shift_{batch}.csv"
        shifts[shifts["cell_id"].isin([f"cell_{c}" for c in cells])].to_csv(
            path, index=False
        )
        paths.append(path)
    return paths


def run(run_component, shift_paths, tmp_path, *extra, min_coverage=10):
    ranked, top = tmp_path / "ranked.csv", tmp_path / "top.csv"
    run_component(
        [
            "--input",
            ";".join(str(p) for p in shift_paths),
            "--output_ranked",
            ranked,
            "--output_top",
            top,
            "--min_coverage",
            str(min_coverage),
            "--n_random",
            "100",
            *extra,
        ]
    )
    return pd.read_csv(ranked), pd.read_csv(top)


def test_ranking_and_top_table(run_component, shift_paths, tmp_path):
    ranked, top = run(run_component, shift_paths, tmp_path)
    assert list(ranked.columns) == [
        "rank",
        "gene_id",
        "gene_name",
        "median_shift",
        "pvalue",
        "n_cells",
    ]
    # RARE has the largest shift but fails the coverage filter
    assert list(ranked["gene_name"]) == ["UP", "ZERO", "DOWN"]
    assert list(ranked["rank"]) == [1, 2, 3]
    assert list(ranked["n_cells"]) == [N_CELLS] * 3
    assert ranked.loc[0, "pvalue"] < 1e-6
    assert ranked.loc[2, "pvalue"] > 0.5
    # only UP is significant from rank 1 onwards
    assert list(top["gene_name"]) == ["UP"]


def test_top_n_and_disease_shift(run_component, shift_paths, tmp_path):
    ranked, top = run(
        run_component,
        shift_paths,
        tmp_path,
        "--shift_column",
        "shift_disease",
        "--top_n",
        "2",
    )
    assert list(ranked["gene_name"]) == ["DOWN", "ZERO", "UP"]
    assert list(top["gene_name"]) == ["DOWN", "ZERO"]


def test_nothing_covered_writes_empty_tables(run_component, shift_paths, tmp_path):
    ranked, top = run(run_component, shift_paths, tmp_path, min_coverage=1000)
    assert ranked.empty and top.empty
    assert "gene_name" in ranked.columns


def test_missing_column_errors(run_component, shift_paths, tmp_path):
    with pytest.raises(CalledProcessError) as err:
        run(run_component, shift_paths, tmp_path, "--shift_column", "shift_unstim")
    assert re.search(r"\['shift_unstim'\] not found", err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
