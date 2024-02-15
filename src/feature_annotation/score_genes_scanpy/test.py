import pytest
import sys
import mudata as mu
import subprocess
import re

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def gene_list_file(tmp_path):
    result = tmp_path / "s_genes.txt"
    gene_list = ["UBE2C", "BIRC5", "TPX2"]
    with result.open('w') as open_gene_list_file:
        open_gene_list_file.write("\n".join(gene_list))
    return result


def test_cell_scoring(run_component, tmp_path):

    output_file = tmp_path / "output.h5mu"

    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--input_layer", "log_normalized",
        "--var_gene_names", "gene_symbol",
        "--gene_list", "UBE2C",
        "--gene_list", "BIRC5",
        "--gene_list", "TPX2",
        "--output", output_file,
        "--obs_score", 'cell_cycle_score'
    ])

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = ["cell_cycle_score"]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, \
            f"could not find columns .mod['rna'].obs['{col}']"


def test_cell_scoring_with_alternative_args(run_component, tmp_path, gene_list_file):
    output_file = tmp_path / "output_new.h5mu"

    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--input_layer", "log_normalized",
        "--var_gene_names", "gene_symbol",
        "--gene_list_file", gene_list_file,
        "--output", output_file,
        "--obs_score", 'cell_cycle_score'
    ])

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = ["cell_cycle_score"]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, \
            f"could not find columns mdata.mod['rna'].obs['{col}']"


def test_fail(run_component, tmp_path):
    output_file = tmp_path / "output_newest.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as e_info:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--input_layer", "log_normalized",
            "--var_gene_names", "gene_symbol",
            "--gene_list", "a_gene_name_that_does_not_exist",
            "--output", output_file
        ])

        assert e_info.value.returncode != 0
        expected_error = r"The following genes are missing from the input dataset: [a_gene_name_that_does_not_exist]"
        assert re.search(expected_error, e_info.value.stdout.decode('utf-8')) is not None, \
            f"expected error message not found in {e_info.value.stdout.decode('utf-8')}"

    assert not output_file.exists(), f"output file should not exist: {output_file}"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
