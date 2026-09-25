import sys
import subprocess
import pytest
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "target/executable/interpret/gseapy/gseapy",
    "resources_dir": "src/interpret/gseapy/",
}
## VIASH END


@pytest.fixture
def test_data(tmp_path):
    """Create a synthetic DE table and two local GMT files in a temp directory."""
    rng = np.random.default_rng(0)
    n_genes = 150
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]

    lfc = rng.normal(0, 2, size=n_genes)
    pvals = rng.uniform(0, 1, size=n_genes)
    # simple BH correction
    rank = np.argsort(pvals).argsort() + 1
    padj = np.clip(pvals * n_genes / rank, 0, 1)
    de = pd.DataFrame(
        {
            "log2FoldChange": lfc,
            "pvalue": pvals,
            "padj": padj,
        },
        index=pd.Index(gene_names, name="gene"),
    )
    csv_path = tmp_path / "deseq2_results.csv"
    de.to_csv(str(csv_path))

    # tiny local GMTs (no internet needed)
    gmt_path = tmp_path / "test_gene_sets.gmt"
    with open(str(gmt_path), "w") as fh:
        fh.write("SET_A\tna\t" + "\t".join(gene_names[:20]) + "\n")
        fh.write("SET_B\tna\t" + "\t".join(gene_names[20:50]) + "\n")
        fh.write("SET_C\tna\t" + "\t".join(gene_names[50:80]) + "\n")

    gmt2_path = tmp_path / "other_gene_sets.gmt"
    with open(str(gmt2_path), "w") as fh:
        fh.write("SET_D\tna\t" + "\t".join(gene_names[80:110]) + "\n")
        fh.write("SET_E\tna\t" + "\t".join(gene_names[110:140]) + "\n")

    return {
        "csv": str(csv_path),
        "gmt": str(gmt_path),
        "gmt2": str(gmt2_path),
        "gene_names": gene_names,
    }


def test_prerank(run_component, tmp_path, test_data):
    output = tmp_path / "enrichment.csv"

    run_component(
        [
            "--input",
            test_data["csv"],
            "--method",
            "prerank",
            "--gene_sets_file",
            test_data["gmt"],
            "--permutation_num",
            "10",
            "--min_size",
            "5",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "Output CSV not created"
    res = pd.read_csv(str(output))
    assert len(res) > 0, "No enrichment terms returned"
    assert list(res.columns[:2]) == ["gene_set_library", "method"], (
        f"Unexpected leading columns: {list(res.columns[:2])}"
    )
    assert set(res["gene_set_library"]) == {"test_gene_sets"}
    assert set(res["method"]) == {"prerank"}
    assert "Term" in res.columns, f"Missing 'Term' column: {list(res.columns)}"


def test_two_libraries_one_file(run_component, tmp_path, test_data):
    """Several gene set libraries end up in one table, not one file each."""
    output = tmp_path / "enrichment_two.csv"

    run_component(
        [
            "--input",
            test_data["csv"],
            "--method",
            "prerank",
            "--gene_sets_file",
            test_data["gmt"],
            "--gene_sets_file",
            test_data["gmt2"],
            "--permutation_num",
            "10",
            "--min_size",
            "5",
            "--output",
            str(output),
        ]
    )

    res = pd.read_csv(str(output))
    assert set(res["gene_set_library"]) == {"test_gene_sets", "other_gene_sets"}, (
        f"Both libraries should be present: {sorted(set(res['gene_set_library']))}"
    )
    # and the rows of each library are the ones that library produced
    assert res.groupby("gene_set_library").size().min() > 0


def test_ora(run_component, tmp_path, test_data):
    output = tmp_path / "enrichment_ora.csv"

    run_component(
        [
            "--input",
            test_data["csv"],
            "--method",
            "ora",
            "--gene_sets_file",
            test_data["gmt"],
            "--pval_threshold",
            "1.0",  # accept all genes so the list is never empty
            "--fc_threshold",
            "0.0",
            "--min_size",
            "5",
            "--output",
            str(output),
        ]
    )

    res = pd.read_csv(str(output))
    assert set(res["method"]) == {"ora"}
    assert len(res) > 0


def test_ora_no_significant_genes_warns(run_component, tmp_path, test_data):
    """No significant genes warns and writes an empty table instead of failing."""
    output = tmp_path / "enrichment_empty.csv"

    stdout = run_component(
        [
            "--input",
            test_data["csv"],
            "--method",
            "ora",
            "--gene_sets_file",
            test_data["gmt"],
            "--pval_threshold",
            "0.0",
            "--min_size",
            "5",
            "--output",
            str(output),
        ]
    )

    assert "No significant genes" in stdout.decode("utf-8")
    res = pd.read_csv(str(output))
    assert len(res) == 0
    assert list(res.columns[:2]) == ["gene_set_library", "method"]


def test_missing_fc_column_fails(run_component, tmp_path, test_data):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                test_data["csv"],
                "--method",
                "prerank",
                "--gene_sets_file",
                test_data["gmt"],
                "--fc_column",
                "not_a_column",
                "--permutation_num",
                "10",
                "--min_size",
                "5",
                "--output",
                str(tmp_path / "enrichment_missing_col.csv"),
            ]
        )
    assert "--fc_column 'not_a_column' not found" in err.value.stdout.decode("utf-8")


def test_non_numeric_fc_column_fails(run_component, tmp_path, test_data):
    """A text column where a fold change is expected must be rejected, not ranked."""
    de = pd.read_csv(test_data["csv"], index_col=0)
    de["log2FoldChange"] = "not_a_number"
    csv_path = tmp_path / "deseq2_text.csv"
    de.to_csv(str(csv_path))

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(csv_path),
                "--method",
                "prerank",
                "--gene_sets_file",
                test_data["gmt"],
                "--permutation_num",
                "10",
                "--min_size",
                "5",
                "--output",
                str(tmp_path / "enrichment_text_col.csv"),
            ]
        )
    assert "is not numeric" in err.value.stdout.decode("utf-8")


def test_no_gene_sets_fails(run_component, tmp_path, test_data):
    """Neither --gene_sets nor --gene_sets_file given."""
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                test_data["csv"],
                "--method",
                "prerank",
                "--permutation_num",
                "10",
                "--min_size",
                "5",
                "--output",
                str(tmp_path / "enrichment_no_gene_sets.csv"),
            ]
        )
    assert "No gene sets provided" in err.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
