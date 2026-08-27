import sys
import subprocess
import pytest
import mudata
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "./target/executable/interpret/cellphonedb_search/",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

# Minimal, hand-crafted CellphoneDB result tables matching the exact schema
# search_analysis_results expects, covering 3 interactions across 3 cell type
# pairs so each query dimension (gene, cell type, classification, score) can
# be exercised independently.
SIGNIFICANT_MEANS = pd.DataFrame(
    {
        "id_cp_interaction": ["CPI-001", "CPI-002", "CPI-003"],
        "interacting_pair": ["GENE1_GENE2", "GENE3_GENE4", "GENE5_GENE6"],
        "partner_a": ["simple:GENE1", "simple:GENE3", "simple:GENE5"],
        "partner_b": ["simple:GENE2", "simple:GENE4", "simple:GENE6"],
        "gene_a": ["GENE1", "GENE3", "GENE5"],
        "gene_b": ["GENE2", "GENE4", "GENE6"],
        "directionality": ["Ligand-Receptor"] * 3,
        "classification": ["Cytokine", "Cytokine", "Chemokine"],
        "Tcell|Bcell": [0.5, np.nan, np.nan],
        "Bcell|Tcell": [np.nan, 0.8, np.nan],
        "Myeloid|Myeloid": [np.nan, np.nan, 1.2],
    }
)

DECONVOLUTED = pd.DataFrame(
    {
        "gene_name": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6"],
        "id_cp_interaction": [
            "CPI-001",
            "CPI-001",
            "CPI-002",
            "CPI-002",
            "CPI-003",
            "CPI-003",
        ],
    }
)

INTERACTION_SCORES = pd.DataFrame(
    {
        "id_cp_interaction": ["CPI-001", "CPI-002", "CPI-003"],
        "Tcell|Bcell": [50, 10, 5],
        "Bcell|Tcell": [10, 80, 5],
        "Myeloid|Myeloid": [5, 5, 90],
    }
)


def prepare_input(tmp_path, with_interaction_scores=True):
    data = mudata.read_h5mu(input_file)
    data.mod["rna"].uns["cellphonedb_significant_means"] = SIGNIFICANT_MEANS
    data.mod["rna"].uns["cellphonedb_deconvoluted"] = DECONVOLUTED
    if with_interaction_scores:
        data.mod["rna"].uns["cellphonedb_interaction_scores"] = INTERACTION_SCORES
    prepared_path = tmp_path / "input_with_cellphonedb_results.h5mu"
    data.write_h5mu(str(prepared_path))
    return prepared_path


def test_search_by_gene(run_component, tmp_path):
    input_path = prepare_input(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--modality",
            "rna",
            "--query_genes",
            "GENE1",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    result = output_data.mod["rna"].uns["cellphonedb_search_results"]
    assert list(result["interacting_pair"]) == ["GENE1_GENE2"]


def test_search_by_cell_types(run_component, tmp_path):
    # Cell type filters only narrow an already gene/interaction/classification/
    # score-selected set of interactions down to matching cell type pair columns
    # - they do not select interactions on their own - so this combines a broad
    # gene query (matching all 3 interactions) with a narrow cell type pair to
    # verify CPI-003 (Myeloid|Myeloid only) gets filtered out.
    input_path = prepare_input(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--modality",
            "rna",
            "--query_genes",
            "GENE1",
            "--query_genes",
            "GENE3",
            "--query_genes",
            "GENE5",
            "--query_cell_types_1",
            "Tcell",
            "--query_cell_types_2",
            "Bcell",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    result = output_data.mod["rna"].uns["cellphonedb_search_results"]
    assert set(result["interacting_pair"]) == {"GENE1_GENE2", "GENE3_GENE4"}
    assert "Myeloid|Myeloid" not in result.columns


def test_search_by_minimum_score(run_component, tmp_path):
    input_path = prepare_input(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--modality",
            "rna",
            "--query_minimum_score",
            "85",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    result = output_data.mod["rna"].uns["cellphonedb_search_results"]
    assert list(result["interacting_pair"]) == ["GENE5_GENE6"]


def test_long_format(run_component, tmp_path):
    input_path = prepare_input(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--modality",
            "rna",
            "--query_genes",
            "GENE1",
            "--query_genes",
            "GENE3",
            "--long_format",
            "--output",
            str(output_path),
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    result = output_data.mod["rna"].uns["cellphonedb_search_results"]
    assert "significant_mean" in result.columns
    assert "interacting_cells" in result.columns


def test_minimum_score_without_interaction_scores_fails(run_component, tmp_path):
    input_path = prepare_input(tmp_path, with_interaction_scores=False)
    output_path = tmp_path / "output.h5mu"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                str(input_path),
                "--modality",
                "rna",
                "--query_minimum_score",
                "50",
                "--output",
                str(output_path),
            ]
        )


def test_missing_cellphonedb_results_fails(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input_file,
                "--modality",
                "rna",
                "--query_genes",
                "GENE1",
                "--output",
                str(output_path),
            ]
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
