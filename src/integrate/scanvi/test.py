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

    expected_obsm_output = "X_scanvi_integrated"
    assert (
        expected_obsm_output in output_rna.obsm
    ), f".obsm['{expected_obsm_output}'] not added\noutput_data: {output_data}"

    # assert that nothing else has changed
    del output_rna.obsm[expected_obsm_output]
    print(output_rna.obs.keys())
    print(input_rna.obs.keys())
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
