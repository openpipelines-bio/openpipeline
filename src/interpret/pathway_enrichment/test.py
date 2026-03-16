import sys
import os
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/interpret/pathway_enrichment/pathway_enrichment",
    "resources_dir": "src/interpret/pathway_enrichment/",
}
## VIASH END


@pytest.fixture
def test_data(tmp_path):
    """Create minimal synthetic h5mu and DESeq2 CSV in a temp directory."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 80, 150
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]

    counts = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype("float32")
    adata = ad.AnnData(X=counts)
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    mdata = mu.MuData({"rna": adata})

    h5mu_path = tmp_path / "input.h5mu"
    mdata.write_h5mu(str(h5mu_path))

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

    # tiny local GMT (3 sets of 20 genes each, no internet needed)
    gmt_path = tmp_path / "test_gene_sets.gmt"
    with open(str(gmt_path), "w") as fh:
        fh.write("SET_A\tna\t" + "\t".join(gene_names[:20]) + "\n")
        fh.write("SET_B\tna\t" + "\t".join(gene_names[20:50]) + "\n")
        fh.write("SET_C\tna\t" + "\t".join(gene_names[50:80]) + "\n")

    return {
        "h5mu": str(h5mu_path),
        "csv": str(csv_path),
        "gmt": str(gmt_path),
        "gene_names": gene_names,
    }


def test_prerank(run_component, tmp_path, test_data):
    output_h5mu = tmp_path / "output.h5mu"
    output_csv_dir = str(tmp_path / "results_prerank")

    run_component(
        [
            "--input", test_data["h5mu"],
            "--input_degenes", test_data["csv"],
            "--method", "prerank",
            "--gene_sets", test_data["gmt"],
            "--permutation_num", "10",
            "--min_size", "5",
            "--output", str(output_h5mu),
            "--output_csv_dir", output_csv_dir,
        ]
    )

    assert output_h5mu.is_file(), "Output h5mu not created"

    mdata = mu.read_h5mu(str(output_h5mu))
    uns = mdata.mod["rna"].uns
    assert "pathway_enrichment" in uns, "uns['pathway_enrichment'] not found"
    assert "test_gene_sets" in uns["pathway_enrichment"], \
        "Gene set key missing in uns['pathway_enrichment']"

    # at least one result term expected
    result_terms = uns["pathway_enrichment"]["test_gene_sets"]
    assert len(next(iter(result_terms.values()))) > 0, "No enrichment terms returned"

    # check CSV was written
    csv_files = os.listdir(output_csv_dir)
    assert any(".prerank.csv" in f for f in csv_files), "No prerank CSV output found"


def test_ora(run_component, tmp_path, test_data):
    output_h5mu = tmp_path / "output.h5mu"
    output_csv_dir = str(tmp_path / "results_ora")

    run_component(
        [
            "--input", test_data["h5mu"],
            "--input_degenes", test_data["csv"],
            "--method", "ora",
            "--gene_sets", test_data["gmt"],
            "--pval_threshold", "1.0",  # accept all genes so the list is never empty
            "--fc_threshold", "0.0",
            "--min_size", "5",
            "--output", str(output_h5mu),
            "--output_csv_dir", output_csv_dir,
        ]
    )

    assert output_h5mu.is_file(), "Output h5mu not created"

    mdata = mu.read_h5mu(str(output_h5mu))
    assert "pathway_enrichment" in mdata.mod["rna"].uns

    csv_files = os.listdir(output_csv_dir)
    assert any(".ora.csv" in f for f in csv_files), "No ORA CSV output found"


def test_uns_key_custom(run_component, tmp_path, test_data):
    output_h5mu = tmp_path / "output.h5mu"

    run_component(
        [
            "--input", test_data["h5mu"],
            "--input_degenes", test_data["csv"],
            "--method", "prerank",
            "--gene_sets", test_data["gmt"],
            "--permutation_num", "10",
            "--min_size", "5",
            "--uns_key", "my_enrichment",
            "--output", str(output_h5mu),
            "--output_csv_dir", str(tmp_path / "results_custom"),
        ]
    )

    mdata = mu.read_h5mu(str(output_h5mu))
    assert "my_enrichment" in mdata.mod["rna"].uns, \
        "Custom uns_key 'my_enrichment' not found"


def test_input_preserved(run_component, tmp_path, test_data):
    """h5mu data (X, var, obs) must be unchanged after adding enrichment results."""
    output_h5mu = tmp_path / "output.h5mu"

    run_component(
        [
            "--input", test_data["h5mu"],
            "--input_degenes", test_data["csv"],
            "--method", "prerank",
            "--gene_sets", test_data["gmt"],
            "--permutation_num", "10",
            "--min_size", "5",
            "--output", str(output_h5mu),
            "--output_csv_dir", str(tmp_path / "results_preserved"),
        ]
    )

    orig = mu.read_h5mu(test_data["h5mu"])
    out = mu.read_h5mu(str(output_h5mu))

    np.testing.assert_array_equal(
        orig.mod["rna"].var_names, out.mod["rna"].var_names
    )
    np.testing.assert_array_equal(
        orig.mod["rna"].obs_names, out.mod["rna"].obs_names
    )
    np.testing.assert_array_equal(
        orig.mod["rna"].X.toarray() if hasattr(orig.mod["rna"].X, "toarray")
        else orig.mod["rna"].X,
        out.mod["rna"].X.toarray() if hasattr(out.mod["rna"].X, "toarray")
        else out.mod["rna"].X,
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
