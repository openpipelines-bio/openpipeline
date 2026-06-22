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
    assert output_rna.n_obs == input_rna.n_obs, (
        f"Number of observations changed\noutput_data: {output_data}"
    )
    assert output_rna.n_vars == input_rna.n_vars, (
        f"Number of variables changed\noutput_data: {output_data}"
    )

    assert "X_scanvi_integrated" in output_rna.obsm, (
        f".obsm['X_scanvi_integrated'] not added\noutput_data: {output_data}"
    )

    assert "scanvi_pred" in output_rna.obs, (
        f".obs['scanvi_pred'] not added\noutput_data: {output_data}"
    )

    assert "scanvi_proba" in output_rna.obs, (
        f".obs['scanvi_proba'] not added\noutput_data: {output_data}"
    )

    predictions = output_data.mod["rna"].obs["scanvi_pred"]
    probabilities = output_data.mod["rna"].obs["scanvi_proba"]

    assert predictions.dtype == "category", (
        "Calculated predictions should be category dtype"
    )
    assert not all(predictions.isna()), "Not all predictions should be NA"
    assert probabilities.dtype == "float32", (
        "Calculated probabilities should be float32 dtype"
    )
    assert all(0 <= value <= 1 for value in probabilities), (
        ".obs at celltypist_probability has values outside the range [0, 1]"
    )

    # assert that nothing else has changed
    del output_rna.obsm["X_scanvi_integrated"]
    del output_rna.obs["scanvi_pred"]
    del output_rna.obs["scanvi_proba"]

    assert_equal(input_rna, output_rna)


def test_scanvi_does_not_predict_unused_categories(
    run_component, random_path, random_h5mu_path, write_mudata_to_file
):
    # When the labels column is a Categorical that declares categories with zero
    # observations (e.g. when a reference is subset from a larger atlas, since
    # subsetting rows does not prune the categorical dtype), scANVI must not
    # predict those absent labels. Unused categories otherwise add untrained
    # output units to the classifier head that can win the inference argmax and
    # emit predictions for labels no reference cell carries.
    import scvi

    input_mdata = mudata.read_h5mu(input_file)
    adata = input_mdata.mod["rna"]

    labels = adata.obs["cell_ontology_class"].astype("category")
    present_labels = set(labels.astype(str).unique())

    phantom_labels = [f"phantom_label_{i}" for i in range(10)]
    adata.obs["cell_ontology_class"] = labels.cat.add_categories(phantom_labels)
    assert not present_labels.intersection(phantom_labels), (
        "Phantom labels must not already be present in the reference"
    )

    input_file_with_phantom = write_mudata_to_file(input_mdata)
    output_file = random_h5mu_path()
    output_model = random_path()

    args = [
        "--input",
        str(input_file_with_phantom),
        "--modality",
        "rna",
        "--obs_labels",
        "cell_ontology_class",
        "--scvi_model",
        scvi_model,
        "--output",
        str(output_file),
        "--output_model",
        str(output_model),
        "--max_epochs",
        "5",
        "--output_compression",
        "gzip",
    ]

    run_component(args)

    # Phantom labels must not be present in trained model architecture.
    registry = scvi.model.SCANVI.load_registry(str(output_model))
    registered_labels = set(
        registry["field_registries"]["labels"]["state_registry"]["categorical_mapping"]
    )
    assert not registered_labels.intersection(phantom_labels), (
        "scANVI registered labels absent from the reference: "
        f"{registered_labels.intersection(phantom_labels)}"
    )

    # Phantom labels should not be present in output data.
    output_data = mudata.read_h5mu(output_file)
    predictions = set(output_data.mod["rna"].obs["scanvi_pred"].astype(str).unique())
    assert predictions.issubset(present_labels), (
        "scANVI predicted labels absent from the reference: "
        f"{predictions - present_labels}"
    )


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
