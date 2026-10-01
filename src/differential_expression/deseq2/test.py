import json
import subprocess
import sys
import pytest
import mudata as mu
import numpy as np
import pandas as pd
import re


## VIASH START
meta = {"resources_dir": "resources_test/"}
## VIASH END

sys.path.append(meta["resources_dir"])


@pytest.fixture
def pseudobulk_test_data_path():
    """Path to the pseudobulk test data"""
    return f"{meta['resources_dir']}/TS_Blood_filtered_pseudobulk.h5mu"


def test_simple_deseq2_execution(run_component, tmp_path, pseudobulk_test_data_path):
    """Test basic DESeq2 execution with minimal parameters"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"

    # Assert one CSV file was found
    csv_files = list(output_dir.glob("deseq2_analysis*.csv"))  # Default prefix
    assert len(csv_files) == 1, "Should find exactly one CSV file in output folder"
    output_file = output_dir / "deseq2_analysis.csv"  # Default prefix
    assert output_file.exists(), "Output CSV file does not exist"

    # Check the output file
    results = pd.read_csv(output_file)
    expected_columns = [
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "significant",
        "gene_id",
        "contrast",
        "comparison_group",
        "control_group",
        "abs_log2FoldChange",
    ]
    assert all(col in results.columns for col in expected_columns), (
        f"Expected columns {expected_columns} not found"
    )

    expected_float_cols = [
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "abs_log2FoldChange",
    ]
    float_cols = results.select_dtypes(include=["float"]).columns.tolist()
    assert all(col in float_cols for col in expected_float_cols), (
        f"Expected float columns {expected_float_cols} not found"
    )

    expected_obj_cols = ["gene_id", "contrast", "comparison_group", "control_group"]
    obj_cols = results.select_dtypes(include=["object"]).columns.tolist()
    assert all(col in obj_cols for col in expected_obj_cols), (
        f"Expected object columns {expected_obj_cols} not found"
    )


def test_simple_deseq2_with_cell_group(
    run_component, tmp_path, pseudobulk_test_data_path
):
    """Test DESeq2 execution with cell groups - should create separate files"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--obs_cell_group",
            "cell_type",
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"

    # Check that multiple CSV files exist (one per cell group)
    csv_files = list(output_dir.glob("deseq2_analysis_*.csv"))  # Default prefix
    assert len(csv_files) > 0, "No cell group-specific CSV files found"

    # Check each file has expected structure
    for csv_file in csv_files:
        results = pd.read_csv(csv_file)
        expected_columns = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "significant",
            "gene_id",
            "contrast",
            "comparison_group",
            "control_group",
            "abs_log2FoldChange",
        ]
        assert all(col in results.columns for col in expected_columns), (
            f"Expected columns {expected_columns} not found in {csv_file}"
        )
        assert len(results) > 0, f"No results found in {csv_file}"

        # Check that all rows have the same cell_type value
        assert results["cell_type"].nunique() == 1, (
            f"Multiple cell types found in {csv_file}"
        )

        expected_float_cols = [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "abs_log2FoldChange",
        ]
        float_cols = results.select_dtypes(include=["float"]).columns.tolist()
        assert all(col in float_cols for col in expected_float_cols), (
            f"Expected float columns {expected_float_cols} not found"
        )

        expected_obj_cols = ["gene_id", "contrast", "comparison_group", "control_group"]
        obj_cols = results.select_dtypes(include=["object"]).columns.tolist()
        assert all(col in obj_cols for col in expected_obj_cols), (
            f"Expected object columns {expected_obj_cols} not found"
        )


def test_complex_design_formula(run_component, tmp_path, pseudobulk_test_data_path):
    """Test DESeq2 with complex design formula accounting for multiple factors"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--design_formula",
            "~ disease + treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
            "--p_adj_threshold",
            "0.1",
            "--log2fc_threshold",
            "0.5",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"
    output_file = output_dir / "deseq2_analysis.csv"  # Default prefix
    assert output_file.exists(), "Output CSV file does not exist"

    results = pd.read_csv(output_file)
    assert "significant" in results.columns, "Significance column not found"
    assert len(results) > 0, "No results found for complex design"


def test_complex_design_formula_with_cell_groups(
    run_component, tmp_path, pseudobulk_test_data_path
):
    """Test that without cell group specified, a single CSV is created in output directory"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--design_formula",
            "~ treatment + disease",
            "--obs_cell_group",
            "cell_type",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"

    # Check that cell group specific files exist
    csv_files = list(output_dir.glob("deseq2_analysis_*.csv"))  # Default prefix
    assert len(csv_files) >= 1, "Could not find cell group-specific files"

    # Check the main file structure
    for csv_file in csv_files:
        results = pd.read_csv(csv_file)
        expected_columns = ["log2FoldChange", "pvalue", "padj", "significant"]
        assert all(col in results.columns for col in expected_columns), (
            f"Expected columns {expected_columns} not found in {csv_file}"
        )
        assert len(results) > 0, f"No results found in {csv_file}"
        assert results["cell_type"].nunique() == 1, (
            f"Multiple cell types found in {csv_file} - should be one per file"
        )


def test_invalid_contrast_column(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that invalid contrast column raises appropriate error"""
    output_dir = tmp_path / "deseq2_output"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                pseudobulk_test_data_path,
                "--output_dir",
                str(output_dir),
                "--design_formula",
                "~ treatment",
                "--contrast_column",
                "nonexistent_column",
                "--contrast_values",
                "group1",
                "--contrast_values",
                "group2",
            ]
        )

    assert re.search(
        r"Missing required columns in metadata: nonexistent_column",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


def test_invalid_design_column(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that invalid contrast column raises appropriate error"""
    output_dir = tmp_path / "deseq2_output"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                pseudobulk_test_data_path,
                "--output_dir",
                str(output_dir),
                "--design_formula",
                "malformed formula",
                "--contrast_column",
                "treatment",
                "--contrast_values",
                "ctrl",
                "--contrast_values",
                "stim",
            ]
        )

    assert re.search(
        r"Invalid design formula: 'malformed formula'",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


def test_custom_output_prefix(run_component, tmp_path, pseudobulk_test_data_path):
    """Test custom output prefix functionality"""
    output_dir = tmp_path / "deseq2_output"
    custom_prefix = "my_custom_analysis"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--output_prefix",
            custom_prefix,
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"

    # Should have exactly one CSV file with custom prefix
    output_file = output_dir / f"{custom_prefix}.csv"
    assert output_file.exists(), (
        f"Custom prefix output CSV file {output_file} does not exist"
    )

    # Check the main file structure
    results = pd.read_csv(output_file)
    expected_columns = ["log2FoldChange", "pvalue", "padj", "significant"]
    assert all(col in results.columns for col in expected_columns), (
        f"Expected columns {expected_columns} not found"
    )
    assert len(results) > 0, "No results found in output"


def test_custom_output_prefix_with_cell_groups(
    run_component, tmp_path, pseudobulk_test_data_path
):
    """Test custom output prefix with cell groups"""
    output_dir = tmp_path / "deseq2_output"
    custom_prefix = "celltype_analysis"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--output_prefix",
            custom_prefix,
            "--obs_cell_group",
            "cell_type",
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert output_dir.exists(), "Output directory does not exist"

    # Check that multiple CSV files exist with custom prefix
    csv_files = list(output_dir.glob(f"{custom_prefix}_*.csv"))
    assert len(csv_files) > 0, (
        f"No cell group-specific CSV files found with prefix {custom_prefix}"
    )

    # Check each file has expected structure
    for csv_file in csv_files:
        results = pd.read_csv(csv_file)
        expected_columns = [
            "log2FoldChange",
            "pvalue",
            "padj",
            "significant",
            "cell_type",
        ]
        assert all(col in results.columns for col in expected_columns), (
            f"Expected columns {expected_columns} not found in {csv_file}"
        )
        assert len(results) > 0, f"No results found in {csv_file}"

        # Check that all rows have the same cell_type value
        assert results["cell_type"].nunique() == 1, (
            f"Multiple cell types found in {csv_file}"
        )


def test_export_normalized_counts(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that --export_normalized_counts writes the sample table, count tables and metadata"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--design_formula",
            "~ disease + treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
            "--export_normalized_counts",
        ]
    )

    expected_files = {
        "deseq2_analysis.csv",
        "deseq2_analysis_samples.csv",
        "deseq2_analysis_normalized_counts.csv",
        "deseq2_analysis_vst.csv",
        "deseq2_analysis_metadata.json",
    }
    found_files = {f.name for f in output_dir.iterdir()}
    assert found_files == expected_files, (
        f"Expected {expected_files}, found {found_files}"
    )

    mod = mu.read_h5mu(pseudobulk_test_data_path)["rna"]
    counts = pd.DataFrame(
        np.round(np.asarray(mod.X)), index=mod.obs_names, columns=mod.var_names
    )

    samples = pd.read_csv(output_dir / "deseq2_analysis_samples.csv")
    assert list(samples.columns) == [
        "sample",
        "disease",
        "treatment",
        "size_factor",
        "library_size",
    ]
    assert list(samples["sample"]) == list(mod.obs_names)
    assert list(samples["treatment"]) == list(mod.obs["treatment"].astype(str))
    assert np.all(np.isfinite(samples["size_factor"]))
    assert np.all(samples["size_factor"] > 0)
    np.testing.assert_array_equal(samples["library_size"], counts.sum(axis=1))

    normalized = pd.read_csv(output_dir / "deseq2_analysis_normalized_counts.csv")
    vst = pd.read_csv(output_dir / "deseq2_analysis_vst.csv")
    for table in [normalized, vst]:
        assert list(table.columns) == ["gene_id"] + list(mod.obs_names)
        assert list(table["gene_id"]) == list(mod.var_names)

    expected_normalized = counts.T / samples["size_factor"].to_numpy()
    np.testing.assert_allclose(
        normalized.drop(columns="gene_id").to_numpy(),
        expected_normalized.to_numpy(),
        rtol=1e-6,
    )

    vst_values = vst.drop(columns="gene_id").to_numpy()
    assert np.all(np.isfinite(vst_values)), "VST values should be finite"
    # VST is a monotonic transformation of the normalized counts within a sample
    first_sample = normalized.columns[1]
    order = np.argsort(normalized[first_sample].to_numpy())
    assert np.all(np.diff(vst[first_sample].to_numpy()[order]) >= -1e-8)

    with open(output_dir / "deseq2_analysis_metadata.json") as f:
        metadata = json.load(f)
    assert metadata["design_formula"] == "~ disease + treatment"
    assert metadata["contrasts"] == [
        {
            "name": "stim_vs_ctrl",
            "comparison_group": "stim",
            "control_group": "ctrl",
        }
    ]
    assert metadata["cell_group"] is None
    assert metadata["n_samples"] == mod.n_obs
    assert metadata["n_genes"] == mod.n_vars
    assert metadata["variance_stabilization"]["blind"] is True
    assert metadata["versions"]["DESeq2"], "DESeq2 version should be recorded"


def test_export_normalized_counts_with_cell_groups(
    run_component, tmp_path, pseudobulk_test_data_path
):
    """Test that each cell group gets its own sample table, count tables and metadata"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--obs_cell_group",
            "cell_type",
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
            "--export_normalized_counts",
        ]
    )

    mod = mu.read_h5mu(pseudobulk_test_data_path)["rna"]
    results_files = [
        f
        for f in output_dir.glob("deseq2_analysis_*.csv")
        if not f.stem.endswith(("_samples", "_normalized_counts", "_vst"))
    ]
    assert len(results_files) == mod.obs["cell_type"].nunique()

    all_samples = []
    for results_file in results_files:
        stem = results_file.stem
        cell_type = pd.read_csv(results_file)["cell_type"].iloc[0]
        group_samples = list(mod.obs_names[mod.obs["cell_type"] == cell_type])

        samples = pd.read_csv(output_dir / f"{stem}_samples.csv")
        assert list(samples["sample"]) == group_samples
        assert set(samples["cell_type"]) == {cell_type}
        all_samples += group_samples

        for suffix in ["normalized_counts", "vst"]:
            table = pd.read_csv(output_dir / f"{stem}_{suffix}.csv")
            assert list(table.columns) == ["gene_id"] + group_samples
            assert len(table) == mod.n_vars

        with open(output_dir / f"{stem}_metadata.json") as f:
            metadata = json.load(f)
        assert metadata["cell_group"] == {"column": "cell_type", "value": cell_type}
        assert metadata["n_samples"] == len(group_samples)

    assert sorted(all_samples) == sorted(mod.obs_names), (
        "Every sample is in one cell group"
    )


def test_var_gene_symbols(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that --var_gene_symbols adds a gene_name column and keeps gene_id unique"""
    mdata = mu.read_h5mu(pseudobulk_test_data_path)
    # Non-unique symbols, as in real annotations
    symbols = [f"SYMBOL{i % 1000}" for i in range(mdata["rna"].n_vars)]
    mdata["rna"].var["symbol"] = symbols
    input_path = tmp_path / "input_with_symbols.h5mu"
    mdata.write_h5mu(input_path)
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            str(input_path),
            "--output_dir",
            str(output_dir),
            "--var_gene_symbols",
            "symbol",
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
            "--export_normalized_counts",
        ]
    )

    expected_symbols = pd.Series(symbols, index=mdata["rna"].var_names)

    results = pd.read_csv(output_dir / "deseq2_analysis.csv")
    assert results["gene_id"].is_unique
    assert set(results["gene_id"]) == set(mdata["rna"].var_names)
    assert list(results["gene_name"]) == list(expected_symbols[results["gene_id"]])

    for suffix in ["normalized_counts", "vst"]:
        table = pd.read_csv(output_dir / f"deseq2_analysis_{suffix}.csv")
        assert list(table.columns[:2]) == ["gene_id", "gene_name"]
        assert list(table["gene_name"]) == symbols

    with open(output_dir / "deseq2_analysis_metadata.json") as f:
        assert json.load(f)["var_gene_symbols"] == "symbol"


def test_invalid_var_gene_symbols(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that a missing --var_gene_symbols column raises an error"""
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                pseudobulk_test_data_path,
                "--output_dir",
                str(tmp_path / "deseq2_output"),
                "--var_gene_symbols",
                "nonexistent_column",
                "--design_formula",
                "~ treatment",
                "--contrast_column",
                "treatment",
                "--contrast_values",
                "stim",
                "--contrast_values",
                "ctrl",
            ]
        )

    assert re.search(
        r"var_gene_symbols 'nonexistent_column' not found",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


def test_no_export_by_default(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that without --export_normalized_counts only the results CSV is written"""
    output_dir = tmp_path / "deseq2_output"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output_dir",
            str(output_dir),
            "--design_formula",
            "~ treatment",
            "--contrast_column",
            "treatment",
            "--contrast_values",
            "stim",
            "--contrast_values",
            "ctrl",
        ]
    )

    assert [f.name for f in output_dir.iterdir()] == ["deseq2_analysis.csv"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
