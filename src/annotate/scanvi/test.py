import pytest
import re
import subprocess
from pathlib import Path

import mudata
from anndata.tests.helpers import assert_equal

## VIASH START
meta = {
    "executable": "./target/executable/integrate/scanvi/scanvi",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

import sys

sys.path.append(meta["resources_dir"])

input_file = f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"
input_file_2 = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
scvi_model = f"{meta['resources_dir']}/scvi_model"


def test_scanvi(run_component):

    args = [
        "--input",
        input_file,
        "--modality",
        "rna",
        "--obs_labels",
        "cell_ontology_class",
        "--scvi_model",
        scvi_model,
        "--output",
        "output.h5mu",
        "--output_model",
        "scanvi_model",
        "--max_epochs",
        "5",
        "--output_compression",
        "gzip",
    ]

    run_component(args)

    input_rna = mudata.read_h5ad(input_file.strip(), mod="rna")
    
    # check files
    assert Path("output.h5mu").is_file(), "Output file does not exist"
    assert Path("scanvi_model").is_dir()
    assert Path("scanvi_model/model.pt").is_file()

    # check output h5mu
    output_data = mudata.read_h5mu("output.h5mu")
    output_rna = output_data.mod["rna"]
    assert (
        output_rna.n_obs == input_rna.n_obs
    ), f"Number of observations changed\noutput_data: {output_data}"
    assert (
        output_rna.n_vars == input_rna.n_vars
    ), f"Number of variables changed\noutput_data: {output_data}"

    assert (
        "X_scanvi_integrated" in output_rna.obsm
    ), f".obsm['X_scanvi_integrated'] not added\noutput_data: {output_data}"

    assert (
        "scanvi_pred" in output_rna.obs
    ), f".obs['scanvi_pred'] not added\noutput_data: {output_data}"
    
    assert (
        "scanvi_proba" in output_rna.obs
    ), f".obs['scanvi_proba'] not added\noutput_data: {output_data}"
    
    predictions = output_data.mod["rna"].obs["scanvi_pred"]
    probabilities = output_data.mod["rna"].obs["scanvi_proba"]
    
    assert (
        predictions.dtype == "category"
    ), "Calculated predictions should be category dtype"
    assert not all(predictions.isna()), "Not all predictions should be NA"
    assert (
        probabilities.dtype == "float32"
    ), "Calculated probabilities should be float32 dtype"
    assert all(
        0 <= value <= 1 for value in probabilities
    ), ".obs at celltypist_probability has values outside the range [0, 1]"

    # assert that nothing else has changed
    del output_rna.obsm["X_scanvi_integrated"]
    del output_rna.obs["scanvi_pred"]
    del output_rna.obs["scanvi_proba"]

    assert_equal(input_rna, output_rna)


def test_raises_with_noncompatible_input_file(run_component):
    args = [
        "--input",
        input_file_2,
        "--modality",
        "rna",
        "--obs_labels",
        "cell_ontology_class",
        "--scvi_model",
        scvi_model,
        "--output",
        "output.h5mu",
        "--output_model",
        "scanvi_model",
        "--max_epochs",
        "5",
        "--output_compression",
        "gzip",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: Number of vars in `adata_target` not the same as source.",
        err.value.stdout.decode("utf-8"),
    )
    

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
