import os
import subprocess
import mudata as mu
import sys
import pytest
import pandas as pd
import re


## VIASH START
meta = {
    "resources_dir": "resources_test/"
}
## VIASH END

sys.path.append(meta["resources_dir"])


@pytest.fixture
def pseudobulk_test_data_path():
    """Path to the pseudobulk test data"""
    return f"{meta['resources_dir']}/TS_Blood_filtered_pseudobulk.h5mu"


def test_simple_deseq2_execution(run_component, tmp_path, pseudobulk_test_data_path):
    """Test basic DESeq2 execution with minimal parameters"""
    output_path = tmp_path / "deseq2_results.csv"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output",
            str(output_path),
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

    assert output_path.exists(), "Output CSV file does not exist"
    
    # Check the output file
    results = pd.read_csv(output_path)
    expected_columns = ["log2FoldChange", "pvalue", "padj", "significant"]
    assert all(col in results.columns for col in expected_columns), f"Expected columns {expected_columns} not found"
    assert len(results) > 0, "No results found in output"


def test_complex_design_formula(run_component, tmp_path, pseudobulk_test_data_path):
    """Test DESeq2 with complex design formula accounting for multiple factors"""
    output_path = tmp_path / "deseq2_complex_results.csv"

    run_component(
        [
            "--input",
            pseudobulk_test_data_path,
            "--output",
            str(output_path),
            "--design_formula",
            "~ cell_type + disease + treatment",
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

    assert output_path.exists(), "Output CSV file does not exist"
    
    results = pd.read_csv(output_path)
    assert "significant" in results.columns, "Significance column not found"
    assert len(results) > 0, "No results found for complex design"


def test_invalid_contrast_column(run_component, tmp_path, pseudobulk_test_data_path):
    """Test that invalid contrast column raises appropriate error"""
    output_path = tmp_path / "deseq2_invalid_results.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                pseudobulk_test_data_path,
                "--output",
                str(output_path),
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
        r"Missing required columns in metadata: \['nonexistent_column'\]",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))