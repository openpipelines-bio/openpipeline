import sys

import mudata as mu
import pandas as pd
import pytest

## VIASH START
meta = {"resources_dir": "resources_test/"}
## VIASH END

input_path = meta["resources_dir"] + "/TS_Blood_filtered.h5mu"


def test_default_execution(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    markers_path = tmp_path / "markers.csv"

    run_component(
        [
            "--input",
            input_path,
            "--modality",
            "rna",
            "--sel_clustering",
            "cell_type",
            "--output",
            str(output_path),
            "--output_markers",
            str(markers_path),
        ]
    )

    assert output_path.is_file(), "Output MuData file was not created"
    output_mdata = mu.read_h5mu(output_path)
    assert "rna" in output_mdata.mod
    assert "rank_genes_groups" in output_mdata["rna"].uns

    assert markers_path.is_file(), "Output markers CSV was not created"
    results = pd.read_csv(markers_path)
    expected_columns = {
        "group",
        "names",
        "scores",
        "pvals",
        "pvals_adj",
        "logfoldchanges",
    }
    assert expected_columns.issubset(results.columns)
    assert len(results) > 0
    assert "is_marker" not in results.columns


def test_filter_results(run_component, tmp_path):
    output_path = tmp_path / "output_filtered.h5mu"
    markers_path = tmp_path / "markers_filtered.csv"

    run_component(
        [
            "--input",
            input_path,
            "--modality",
            "rna",
            "--sel_clustering",
            "cell_type",
            "--output",
            str(output_path),
            "--output_markers",
            str(markers_path),
            "--filter_results",
            "true",
        ]
    )

    output_mdata = mu.read_h5mu(output_path)
    assert "rank_genes_groups" in output_mdata["rna"].uns

    filtered_key = "rank_genes_groups_filtered"
    assert filtered_key in output_mdata["rna"].uns
    filtered_uns = output_mdata["rna"].uns[filtered_key]
    assert "params" in filtered_uns
    assert "filters" in filtered_uns
    assert "markers" in filtered_uns

    results = pd.read_csv(markers_path)
    assert "is_marker" in results.columns
    assert set(results["is_marker"].unique()).issubset({True, False})


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
