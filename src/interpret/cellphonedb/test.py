import sys
import subprocess
import pytest
import mudata
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "./target/executable/interpret/cellphonedb/",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
groupby = "harmony_integration_leiden_1.0"


def assert_input_unchanged(input_data, output_data):
    np.testing.assert_array_equal(
        output_data.mod["rna"].X.data, input_data.mod["rna"].X.data
    )
    np.testing.assert_array_equal(
        input_data.mod["rna"].var.index, output_data.mod["rna"].var.index
    )


def test_statistical_analysis(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--output_compression",
            "gzip",
            "--modality",
            "rna",
            "--layer",
            "log_normalized",
            "--groupby",
            groupby,
            "--gene_symbol",
            "gene_symbol",
            "--counts_data",
            "hgnc_symbol",
            "--method",
            "statistical_analysis",
            "--iterations",
            "60",
            "--threshold",
            "0.1",
            "--pvalue",
            "0.05",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    input_data = mudata.read_h5mu(input_file)
    output_data = mudata.read_h5mu(output_path)
    assert_input_unchanged(input_data, output_data)

    uns = output_data.mod["rna"].uns
    for key in [
        "cellphonedb_means",
        "cellphonedb_pvalues",
        "cellphonedb_significant_means",
        "cellphonedb_deconvoluted",
        "cellphonedb_deconvoluted_percents",
    ]:
        assert key in uns


def test_analysis(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--layer",
            "log_normalized",
            "--groupby",
            groupby,
            "--gene_symbol",
            "gene_symbol",
            "--counts_data",
            "hgnc_symbol",
            "--method",
            "analysis",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    uns = output_data.mod["rna"].uns
    assert "cellphonedb_means" in uns
    assert "cellphonedb_deconvoluted" in uns
    assert "cellphonedb_deconvoluted_percents" in uns
    assert "cellphonedb_pvalues" not in uns
    assert "cellphonedb_significant_means" not in uns


def test_degs_analysis(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    input_data = mudata.read_h5mu(input_file)
    var = input_data.mod["rna"].var
    obs = input_data.mod["rna"].obs
    clusters = obs[groupby].astype(str).unique()
    degs = pd.DataFrame(
        {
            "cluster": clusters,
            "gene": var["gene_symbol"].astype(str).iloc[: len(clusters)].values,
        }
    )
    degs_file = tmp_path / "degs.tsv"
    degs.to_csv(degs_file, sep="\t", index=False)

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--layer",
            "log_normalized",
            "--groupby",
            groupby,
            "--gene_symbol",
            "gene_symbol",
            "--counts_data",
            "hgnc_symbol",
            "--method",
            "degs_analysis",
            "--degs_file",
            str(degs_file),
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    uns = output_data.mod["rna"].uns
    assert "cellphonedb_relevant_interactions" in uns
    assert "cellphonedb_significant_means" in uns


def test_degs_analysis_without_degs_file_fails(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input_file,
                "--modality",
                "rna",
                "--layer",
                "log_normalized",
                "--groupby",
                groupby,
                "--method",
                "degs_analysis",
                "--output",
                str(output_path),
            ]
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
