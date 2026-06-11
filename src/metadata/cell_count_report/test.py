import pytest
import uuid
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

## VIASH START
meta = {
    "name": "cell_count_report",
    "resources_dir": "resources_test/",
    "executable": "target/executable/metadata/cell_count_report/cell_count_report",
    "config": "src/metadata/cell_count_report/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def write_temp_h5mu(tmp_path):
    def wrapper(test_h5mu):
        test_h5mu_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
        test_h5mu.write_h5mu(test_h5mu_path)
        return test_h5mu_path

    return wrapper


@pytest.fixture
def h5mu():
    # Two samples, 4 cells each. Boolean keep-columns (True = keep) chosen so the
    # per-stage and combined counts are easy to verify by hand.
    obs_names = [f"cell_{i}" for i in range(8)]
    sample_id = ["s1"] * 4 + ["s2"] * 4

    # RNA filter masks (AND of these two columns gives the RNA keep set).
    filter_counts_rna = [True, True, True, False, True, True, False, False]
    filter_quantile_rna = [True, True, False, True, True, False, True, True]
    # rna_keep:                T,    T,    F,     F,    T,    F,     F,     F
    filter_scrublet = [True, False, True, True, True, True, True, False]
    # prot masks
    filter_counts_prot = [True, True, True, True, True, True, True, True]
    filter_quantile_prot = [True, True, True, True, False, True, True, True]

    rna_obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "filter_counts_rna": filter_counts_rna,
            "filter_quantile_rna": filter_quantile_rna,
            "filter_scrublet": filter_scrublet,
        },
        index=obs_names,
    )
    prot_obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "filter_counts_prot": filter_counts_prot,
            "filter_quantile_prot": filter_quantile_prot,
        },
        index=obs_names,
    )
    rna = AnnData(np.ones((8, 3)), obs=rna_obs)
    prot = AnnData(np.ones((8, 2)), obs=prot_obs)
    return MuData({"rna": rna, "prot": prot})


def test_cell_count_report(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "rna",
            "--prot_modality",
            "prot",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--scrublet_filter_column",
            "filter_scrublet",
            "--prot_filter_columns",
            "filter_counts_prot;filter_quantile_prot",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "Output TSV must have been created."
    report = pd.read_csv(output, sep="\t").set_index("sample_id")

    expected_columns = {
        "cell_count_after_rna_filter",
        "cell_count_after_scrublet_filter",
        "cell_count_after_protein_filter",
        "cell_count_after_all_filter",
    }
    assert expected_columns.issubset(set(report.columns))
    assert set(report.index) == {"s1", "s2"}

    # s1: rna_keep = [T,T,F,F] -> 2 ; scrublet = [T,F,T,T] -> 3 ; prot = [T,T,T,T] -> 4
    #     all = rna & scrublet & prot = [T,F,F,F] -> 1
    assert report.loc["s1", "cell_count_after_rna_filter"] == 2
    assert report.loc["s1", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s1", "cell_count_after_protein_filter"] == 4
    assert report.loc["s1", "cell_count_after_all_filter"] == 1

    # s2: rna_keep = [T,F,F,F] -> 1 ; scrublet = [T,T,T,F] -> 3 ; prot = [F,T,T,T] -> 3
    #     all = [F,F,F,F]... cell_4: rna T, scrublet T, prot F -> F ; so all -> 0
    assert report.loc["s2", "cell_count_after_rna_filter"] == 1
    assert report.loc["s2", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s2", "cell_count_after_protein_filter"] == 3
    assert report.loc["s2", "cell_count_after_all_filter"] == 0

    # The combined count never exceeds any single-stage count.
    for sample in report.index:
        for col in (
            "cell_count_after_rna_filter",
            "cell_count_after_scrublet_filter",
            "cell_count_after_protein_filter",
        ):
            assert (
                report.loc[sample, "cell_count_after_all_filter"]
                <= report.loc[sample, col]
            )


def test_cell_count_report_skip_scrublet_and_prot(
    run_component, h5mu, write_temp_h5mu, tmp_path
):
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "rna",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--output",
            output,
        ]
    )
    assert output.is_file()
    report = pd.read_csv(output, sep="\t").set_index("sample_id")
    # No scrublet, no protein modality: those columns are omitted.
    assert "cell_count_after_scrublet_filter" not in report.columns
    assert "cell_count_after_protein_filter" not in report.columns
    assert "cell_count_after_rna_filter" in report.columns
    assert "cell_count_after_all_filter" in report.columns
    # With only the RNA filter, the combined count equals the RNA count.
    assert (
        report["cell_count_after_all_filter"] == report["cell_count_after_rna_filter"]
    ).all()
    assert report.loc["s1", "cell_count_after_rna_filter"] == 2
    assert report.loc["s2", "cell_count_after_rna_filter"] == 1


if __name__ == "__main__":
    exit(pytest.main([__file__]))
