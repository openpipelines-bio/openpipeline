import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_cellranger_multi_to_h5mu/from_cellranger_multi_to_h5mu",
    "resources_dir": "resources_test/",
    "config": "src/convert/from_cellranger_multi_to_h5mu/config.vsh.yaml",
}
## VIASH END

input_anticmv = f"{meta['resources_dir']}/10x_5k_anticmv/processed/10x_5k_anticmv.cellranger_multi.output"
input_lung_crispr = f"{meta['resources_dir']}/10x_5k_lung_crispr/processed/10x_5k_lung_crispr.cellranger_multi.output.output"
input_beam = (
    f"{meta['resources_dir']}/10x_5k_beam/processed/10x_5k_beam.cellranger_multi.output"
)
input_fixed_rna = f"{meta['resources_dir']}/10x_5k_fixed/processed/10x_5k_fixed.cellranger_multi.output"


def test_cellranger_multi_basic(run_component, tmp_path):
    output_dir = tmp_path / "converted"
    output_path_template = output_dir / "*.h5mu"
    samples_csv = tmp_path / "samples.csv"
    # run component
    run_component(
        [
            "--input",
            input_anticmv,
            "--output",
            str(output_path_template),
            "--output_compression",
            "gzip",
            "--sample_csv",
            samples_csv,
        ]
    )
    assert output_dir.is_dir()

    # check output
    samples = [item for item in output_dir.iterdir() if item.is_file()]
    assert len(samples) == 1
    output_path = samples[0]
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ["rna", "prot", "vdj_t"]
    assert list(converted_data.uns.keys()) == ["metrics_cellranger"]
    expected_metrics = [
        "Category",
        "Library Type",
        "Grouped By",
        "Group Name",
        "Metric Name",
        "Metric Value",
    ]
    assert (
        converted_data.uns["metrics_cellranger"].columns.to_list() == expected_metrics
    )
    # Check that a metric that is stored as percentage (e.g "85.69%") is correctly represented
    # as a floating point number
    metrics_df_with_index = converted_data.uns["metrics_cellranger"].set_index(
        ["Metric Name", "Library Type", "Category"]
    )
    percentage = metrics_df_with_index.loc[
        ("Confidently mapped reads in cells", "Gene Expression", "Cells"),
        "Metric Value",
    ]
    assert percentage.iloc[0] == "0.8569"

    thousand_delimited_number = metrics_df_with_index.loc[
        ("Cells", "Gene Expression", "Cells"), "Metric Value"
    ]
    thousand_delimited_number == "3,798"

    smaller_number = metrics_df_with_index.loc[
        ("Median genes per cell", "Gene Expression", "Cells"), "Metric Value"
    ]
    smaller_number == "6"


def test_cellranger_multi_to_h5mu_crispr(run_component, tmp_path):
    output_dir = tmp_path / "converted"
    output_path_template = output_dir / "*.h5mu"
    samples_csv = tmp_path / "samples.csv"

    # run component
    run_component(
        [
            "--input",
            input_lung_crispr,
            "--output",
            str(output_path_template),
            "--output_compression",
            "gzip",
            "--sample_csv",
            samples_csv,
        ]
    )
    assert output_dir.is_dir()

    # check output
    samples = [item for item in output_dir.iterdir() if item.is_file()]
    assert len(samples) == 1
    output_path = samples[0]
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ["rna", "gdo"]
    assert list(converted_data.uns.keys()) == ["metrics_cellranger"]
    assert "perturbation_efficiencies_by_feature" in converted_data.mod["gdo"].uns
    assert "perturbation_efficiencies_by_target" in converted_data.mod["gdo"].uns
    assert "feature_reference" not in converted_data.mod["rna"].uns
    assert "feature_reference" in converted_data.mod["gdo"].uns


def test_cellranger_multi_to_h5mu_beam(run_component, tmp_path):
    output_dir = tmp_path / "converted"
    output_path_template = output_dir / "*.h5mu"
    samples_csv = tmp_path / "samples.csv"

    # run component
    run_component(
        [
            "--input",
            input_beam,
            "--output",
            str(output_path_template),
            "--output_compression",
            "gzip",
            "--sample_csv",
            samples_csv,
        ]
    )
    assert output_dir.is_dir()

    # check output
    samples = [item for item in output_dir.iterdir() if item.is_file()]
    assert len(samples) == 1
    output_path = samples[0]
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ["rna", "antigen", "vdj_t"]
    assert "antigen_specificity_scores_CMV_B0702" in converted_data["antigen"].obsm
    assert "antigen_specificity_scores_Flu_A0201" in converted_data["antigen"].obsm


def test_cellranger_multi_to_h5mu_fixed_rna(run_component, tmp_path):
    output_dir = tmp_path / "converted"
    output_path_template = output_dir / "*.h5mu"
    samples_csv = tmp_path / "samples.csv"

    # run component
    run_component(
        [
            "--input",
            input_fixed_rna,
            "--output",
            str(output_path_template),
            "--output_compression",
            "gzip",
            "--sample_csv",
            samples_csv,
        ]
    )
    assert output_dir.is_dir()

    # check output
    samples = [item for item in output_dir.iterdir() if item.is_file()]
    sample_names = {item.name.removesuffix(".h5mu") for item in samples}
    assert sample_names == {
        "Colorectal_BC3",
        "Liver_BC1",
        "Ovarian_BC2",
        "Pancreas_BC4",
    }
    for output_path in samples:
        converted_data = read_h5mu(output_path)
        assert list(converted_data.mod.keys()) == ["rna", "prot"]


def test_custom_modality(run_component, tmp_path):
    input_dir = f"{meta['resources_dir']}/10x_5k_anticmv/processed_with_custom/10x_5k_anticmv.cellranger_multi.output"
    output_dir = tmp_path / "converted"
    output_path_template = output_dir / "*.h5mu"
    samples_csv = tmp_path / "samples.csv"
    # run component
    run_component(
        [
            "--input",
            input_dir,
            "--output",
            str(output_path_template),
            "--output_compression",
            "gzip",
            "--sample_csv",
            samples_csv,
        ]
    )
    assert output_dir.is_dir()

    # check output
    samples = [item for item in output_dir.iterdir() if item.is_file()]
    assert len(samples) == 1
    output_path = samples[0]
    converted_data = read_h5mu(output_path)
    assert list(converted_data.mod.keys()) == ["rna", "custom", "vdj_t"]
    assert "feature_reference" not in converted_data.mod["custom"].uns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
