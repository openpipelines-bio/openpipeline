import pytest
import subprocess
from mudata import read_h5mu
import re
import sys

## VIASH START
meta = {
    "executable": "./target/docker/scgpt/cross_check/cross_check",
    "resources_dir": "./resources_test/scgpt/",
    "config": "./src/scgpt/cross_check/config.vsh.yaml",
}
## VIASH END

input_path = meta["resources_dir"] + "/Kim2020_Lung_subset_preprocessed.h5mu"
vocab_path = meta["resources_dir"] + "/vocab.json"


def test_cross_check(run_component, random_path):
    output_path = random_path(extension="h5mu")
    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--modality",
        "rna",
        "--vocab_file",
        vocab_path,
        "--output_compression",
        "gzip",
    ]
    run_component(args)

    output_mudata = read_h5mu(output_path)

    # Check added columns
    assert {"gene_name", "id_in_vocab"}.issubset(
        set(output_mudata.mod["rna"].var.columns)
    ), "Gene columns were not added."
    # Check if genes were filtered
    assert sum(output_mudata.mod["rna"].var["id_in_vocab"]) != len(
        output_mudata.mod["rna"].var["id_in_vocab"]
    ), "Genes were not filtered."

    output_hvg_path = random_path(extension="h5mu")
    args_hvg = [
        "--input",
        input_path,
        "--output",
        output_hvg_path,
        "--modality",
        "rna",
        "--var_input",
        "filter_with_hvg",
        "--vocab_file",
        vocab_path,
        "--output_compression",
        "gzip",
    ]

    run_component(args_hvg)

    output_mudata_hvg = read_h5mu(output_hvg_path)
    # Check if genes were filtered based on HVG
    assert sum(output_mudata_hvg.mod["rna"].var["id_in_vocab"]) != len(
        output_mudata_hvg.mod["rna"].var["id_in_vocab"]
    ), "Genes were not filtered."
    assert sum(output_mudata.mod["rna"].var["id_in_vocab"]) != len(
        output_mudata_hvg.mod["rna"].var["id_in_vocab"]
    ), "Genes were not filtered based on HVG."


def test_cross_check_invalid_gene_layer_raises(run_component, random_path):
    output_path = random_path(extension="h5mu")
    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--vocab_file",
        vocab_path,
        "--input_var_gene_names",
        "dummy_var",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: Gene name column 'dummy_var' not found in .mod\['rna'\]\.obs\.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
