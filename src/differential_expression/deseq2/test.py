import subprocess
import sys
import pytest
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
            "--padj_threshold",
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
