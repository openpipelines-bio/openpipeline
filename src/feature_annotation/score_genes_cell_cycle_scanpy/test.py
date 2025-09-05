import pytest
import sys
import mudata as mu
import subprocess
import re

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def test_cell_scoring_cell_cycle(run_component, tmp_path):
    output_file = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--input_layer",
            "log_normalized",
            "--var_gene_names",
            "gene_symbol",
            "--s_genes",
            "MCM5",
            "--s_genes",
            "PCNA",
            "--s_genes",
            "TYMS",
            "--g2m_genes",
            "UBE2C",
            "--g2m_genes",
            "BIRC5",
            "--g2m_genes",
            "TPX2",
            "--output",
            output_file,
            "--obs_phase",
            "my_phase",
            "--obs_s_score",
            "my_s_score",
            "--obs_g2m_score",
            "my_g2m_score",
        ]
    )

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = [
        "my_phase",
        "my_s_score",
        "my_g2m_score",
    ]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, (
            f"could not find columns mdata.mod['rna'].obs['{col}']"
        )


def test_cell_scoring_cell_cycle_with_alternative_args(run_component, tmp_path):
    output_file = tmp_path / "output_new.h5mu"
    g2m_gene_file = tmp_path / "g2m_genes.txt"
    s_gene_file = tmp_path / "s_genes.txt"

    with open(g2m_gene_file, "w") as f:
        f.write("UBE2C\nBIRC5\nTPX2\n")
    with open(s_gene_file, "w") as f:
        f.write("MCM5\nPCNA\nTYMS\n")

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--input_layer",
            "log_normalized",
            "--var_gene_names",
            "gene_symbol",
            "--s_genes_file",
            s_gene_file,
            "--g2m_genes_file",
            g2m_gene_file,
            "--output",
            output_file,
        ]
    )

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = [
        "phase",
        "S_score",
        "G2M_score",
    ]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, (
            f"could not find columns mdata.mod['rna'].obs['{col}']"
        )


def test_cell_scoring_cell_cycle_with_mixed_args(run_component, tmp_path):
    output_file = tmp_path / "output_new.h5mu"
    g2m_gene_file = tmp_path / "g2m_genes.txt"
    s_gene_file = tmp_path / "s_genes.txt"

    with open(g2m_gene_file, "w") as f:
        f.write("UBE2C\nBIRC5\nTPX2\n")
    with open(s_gene_file, "w") as f:
        f.write("MCM5\nPCNA\nTYMS\n")

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--input_layer",
            "log_normalized",
            "--var_gene_names",
            "gene_symbol",
            "--s_genes_file",
            s_gene_file,
            "--s_genes",
            "FEN1",
            "--g2m_genes_file",
            g2m_gene_file,
            "--g2m_genes",
            "TOP2A",
            "--output",
            output_file,
        ]
    )

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = [
        "phase",
        "S_score",
        "G2M_score",
    ]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, (
            f"could not find columns mdata.mod['rna'].obs['{col}']"
        )


def test_fail(run_component, tmp_path):
    output_file = tmp_path / "output_newest.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as e_info:
        run_component(
            [
                "--input",
                input_file,
                "--modality",
                "rna",
                "--input_layer",
                "log_normalized",
                "--var_gene_names",
                "gene_symbol",
                "--s_genes",
                "a_gene_name_that_does_not_exist",
                "--g2m_genes",
                "MCM5",
                "--output",
                output_file,
            ]
        )

    assert e_info.value.returncode != 0
    expected_error = r"The follow genes are missing from the input dataset: {\'a_gene_name_that_does_not_exist\'}"
    assert re.search(expected_error, e_info.value.stdout.decode("utf-8")) is not None, (
        f"expected error message not found in {e_info.value.stdout.decode('utf-8')}"
    )

    assert not output_file.exists(), f"output file should not exist: {output_file}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
