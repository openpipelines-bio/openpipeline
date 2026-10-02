import json
import re
import subprocess
import sys

import pandas as pd
import pytest

## VIASH START
meta = {"resources_dir": "src/report/deseq2_report"}
## VIASH END

sys.path.append(meta["resources_dir"])

DESEQ2_DIR = f"{meta['resources_dir']}/test_data"
OVERALL_DIR = f"{DESEQ2_DIR}/overall"
PER_CELL_TYPE_DIR = f"{DESEQ2_DIR}/per_cell_type"


def read_report_data(html_path):
    """Extract the data object embedded in the report."""
    html = html_path.read_text()
    match = re.search(r"window\.REPORT_DATA = (\{.*?\});</script>", html, re.S)
    assert match, "Report data not found in the HTML"
    return json.loads(match.group(1).replace("<\\/", "</"))


def test_report_overall(run_component, tmp_path):
    """Test a report on an analysis without significant genes"""
    output = tmp_path / "report.html"
    output_data = tmp_path / "report_data.json"

    run_component(
        [
            "--input",
            OVERALL_DIR,
            "--output",
            str(output),
            "--output_data",
            str(output_data),
        ]
    )

    assert output.exists(), "Report was not created"
    html = output.read_text()
    assert "/*__" not in html, "Not all template placeholders were filled"
    assert "plotly.js (cartesian - minified) v2.35.2" in html, "plotly.js not inlined"
    assert not re.search(r"<script[^>]+src=", html), "Report should not load scripts"

    data = json.loads(output_data.read_text())
    assert data == read_report_data(output), "HTML and JSON data differ"

    samples = pd.read_csv(f"{OVERALL_DIR}/deseq2_analysis_samples.csv")
    results = pd.read_csv(f"{OVERALL_DIR}/deseq2_analysis.csv")
    tested = results.dropna(subset=["padj"])

    assert data["meta"]["title"] == "Differential expression report"
    assert data["kpis"]["samples"] == len(samples)
    assert data["kpis"]["padj_threshold"] == 0.05, "Threshold of the DESeq2 run"
    assert data["kpis"]["de_genes"] == int(results["significant"].sum())
    assert data["provenance"]["contrast"] == {
        "column": "treatment",
        "test": "stim",
        "reference": "ctrl",
    }
    assert data["provenance"]["deseq2"]["design_formula"] == "~ cell_type + treatment"

    # Samples ordered reference group first
    groups = data["pca"]["group"]
    assert groups == sorted(groups, key=lambda g: g != "ctrl")
    assert len(data["pca"]["pc1"]) == len(samples)
    assert len(data["corr"]["z"]) == len(samples)

    assert data["heatmap"]["genes"] == [], "No significant genes to show"
    assert data["genes"]["items"] == [], "No key genes without significant genes"
    assert data["summary"]["strongest"] is None
    assert data["pvalues"]["n"] == len(tested)
    assert sum(data["pvalues"]["counts"]) == len(tested)
    assert len(data["table"]) == min(300, len(tested))
    assert [row["padj"] for row in data["table"]] == sorted(
        row["padj"] for row in data["table"]
    )
    assert {m[0] for m in data["methods"]} >= {
        "Design",
        "Differential expression",
        "Expression scale",
    }


def test_report_cell_group_with_options(run_component, tmp_path):
    """Test a per cell group report with significant genes and the display options"""
    output = tmp_path / "report.html"
    output_data = tmp_path / "report_data.json"
    prefix = "deseq2_analysis_erythrocyte"
    results = pd.read_csv(f"{PER_CELL_TYPE_DIR}/{prefix}.csv")
    significant = results[results["significant"]]
    highlight = significant.sort_values("padj")["gene_id"].iloc[0]

    run_component(
        [
            "--input",
            PER_CELL_TYPE_DIR,
            "--input_prefix",
            prefix,
            "--output",
            str(output),
            "--output_data",
            str(output_data),
            "--title",
            "Erythrocytes",
            "--project",
            "Test project",
            "--group_labels",
            "stim=Stimulated",
            "--group_labels",
            "ctrl=Control",
            "--obs_sample_label",
            "cell_type",
            "--obs_sample_label",
            "treatment",
            "--highlight_genes",
            highlight,
            "--highlight_genes",
            "NOT_A_GENE",
            "--n_heatmap_genes",
            "5",
            "--n_table_genes",
            "20",
        ]
    )

    data = json.loads(output_data.read_text())
    assert data["meta"]["title"] == "Erythrocytes"
    assert data["meta"]["project"] == "Test project"
    assert set(data["pca"]["group"]) == {"Stimulated", "Control"}
    assert all(label.startswith("erythrocyte ") for label in data["pca"]["sample"])
    assert data["kpis"]["de_genes"] == len(significant)
    assert len(data["heatmap"]["genes"]) == min(5, len(significant))
    assert len(data["table"]) == 20
    assert [item["gene_id"] for item in data["genes"]["items"]] == [highlight], (
        "Unknown highlight genes are skipped"
    )
    assert any("erythrocyte" in m[1] for m in data["methods"] if m[0] == "Design"), (
        "The cell group should be named in the methods"
    )


def test_obs_pair_and_thresholds(run_component, tmp_path):
    """Test pairing samples and overriding the significance thresholds"""
    output_data = tmp_path / "report_data.json"

    run_component(
        [
            "--input",
            OVERALL_DIR,
            "--output",
            str(tmp_path / "report.html"),
            "--output_data",
            str(output_data),
            "--obs_pair",
            "cell_type",
            "--pair_label",
            "cell type",
            "--p_adj_threshold",
            "1",
            "--log2fc_threshold",
            "0.5",
        ]
    )

    data = json.loads(output_data.read_text())
    results = pd.read_csv(f"{OVERALL_DIR}/deseq2_analysis.csv").dropna(subset=["padj"])
    expected = ((results["padj"] < 1) & (results["log2FoldChange"].abs() > 0.5)).sum()
    assert data["kpis"]["padj_threshold"] == 1
    assert data["kpis"]["de_genes"] == expected
    samples = pd.read_csv(f"{OVERALL_DIR}/deseq2_analysis_samples.csv")
    assert sorted(data["genes"]["pairs"]) == sorted(samples["cell_type"])
    assert data["genes"]["pair_label"] == "cell type"
    # Cell type dominates this pseudobulk data, so a PC follows it
    assert "follows cell type" in data["pca"]["hint"]


def test_missing_input_files(run_component, tmp_path):
    """Test that a wrong prefix gives a clear error"""
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                OVERALL_DIR,
                "--input_prefix",
                "does_not_exist",
                "--output",
                str(tmp_path / "report.html"),
            ]
        )
    assert re.search(
        r"Missing DESeq2 output files", err.value.stdout.decode("utf-8")
    ), err.value.stdout.decode("utf-8")


def test_invalid_obs_pair(run_component, tmp_path):
    """Test that a pairing column outside the DESeq2 sample table gives a clear error"""
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                OVERALL_DIR,
                "--obs_pair",
                "disease",
                "--output",
                str(tmp_path / "report.html"),
            ]
        )
    assert re.search(
        r"Column 'disease' not found in the DESeq2 sample table",
        err.value.stdout.decode("utf-8"),
    ), err.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
