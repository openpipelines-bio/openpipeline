import pandas as pd
import sys
import pytest

## VIASH START
par = {"input": "test/deseq_analysis_classical_monocyte.csv"}
## VIASH END


def test_run():
    results = pd.read_csv(par["input"])

    # Make sure expected columns are present
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

    # Make sure columns have expected dtypes
    assert results["significant"].dtypes == "bool", (
        "Expected 'significant' column to be of type bool"
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
