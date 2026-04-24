import sys
from subprocess import CalledProcessError

import mudata as mu
import pytest

from openpipeline_testutils.asserters import assert_annotation_objects_equal


## VIASH START
meta = {
    "resources_dir": "./resources_test/",
}
## VIASH END


ts_blood_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
pbmc_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3/"
    "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)


@pytest.mark.parametrize("output_compression", [None, "gzip", "lzf"])
def test_clear_multiple_slots(run_component, random_h5mu_path, output_compression):
    temp_output = random_h5mu_path()

    arguments = [
        "--input",
        ts_blood_file,
        "--modality",
        "rna",
        "--slots",
        "obsm",
        "--slots",
        "varm",
        "--output",
        temp_output,
    ]
    if output_compression:
        arguments.extend(["--output_compression", output_compression])

    run_component(arguments)

    output = mu.read_h5mu(temp_output)
    assert len(output.mod["rna"].obsm) == 0
    assert "log_normalized" in output.mod["rna"].layers
    assert set(output.mod) == {"rna"}

    expected = mu.read_h5mu(ts_blood_file)
    expected.mod["rna"].obsm.clear()
    assert_annotation_objects_equal(expected, output)


@pytest.mark.parametrize(
    ("slot", "source_file"),
    [
        ("varm", ts_blood_file),
        ("obsp", ts_blood_file),
        ("varp", ts_blood_file),
        ("uns", pbmc_file),
    ],
)
def test_clear_individual_supported_slots(
    run_component, random_h5mu_path, slot, source_file
):
    temp_output = random_h5mu_path()
    input_data = mu.read_h5mu(source_file)
    before_len = len(getattr(input_data.mod["rna"], slot))

    run_component(
        [
            "--input",
            source_file,
            "--modality",
            "rna",
            "--slots",
            slot,
            "--output",
            temp_output,
        ]
    )

    output = mu.read_h5mu(temp_output)
    assert len(getattr(output.mod["rna"], slot)) == 0
    assert before_len >= 0


def test_multimodal_pbmc_preserves_untouched_modality(run_component, random_h5mu_path):
    temp_output = random_h5mu_path()

    run_component(
        [
            "--input",
            pbmc_file,
            "--modality",
            "rna",
            "--slots",
            "uns",
            "--output",
            temp_output,
        ]
    )

    output = mu.read_h5mu(temp_output)
    expected = mu.read_h5mu(pbmc_file)
    expected.mod["rna"].uns.clear()
    assert set(output.mod) == {"rna", "prot"}
    assert_annotation_objects_equal(expected, output)


def test_missing_modality_raises(run_component, random_h5mu_path):
    temp_output = random_h5mu_path()

    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                pbmc_file,
                "--modality",
                "missing_modality",
                "--slots",
                "obsm",
                "--output",
                temp_output,
            ]
        )

    assert not temp_output.is_file()
    assert "modality" in err.value.stdout.decode("utf-8").lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
