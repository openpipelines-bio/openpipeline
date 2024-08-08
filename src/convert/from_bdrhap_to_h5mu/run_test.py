from os import path
from mudata import read_h5mu
import numpy as np
import pytest
import sys

## VIASH START
meta = {
    'resources_dir': 'resources_test/',
    'config': './src/convert/from_bdrhap_to_h5mu/config.vsh.yaml',
    'executable': './target/executable/convert/from_bdrhap_to_h5mu/from_bdrhap_to_h5mu',
}
## VIASH END

input = meta["resources_dir"] + "/WTA.bd_rhapsody.output_raw"

def test_run(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input", input,
        "--output", output,
        "--id", "foo",
        "--output_compression", "gzip",
    ]
    run_component(cmd_pars)

    # check if file exists
    assert path.exists(output), "No output was created."

    # read it with scanpy
    data = read_h5mu(output)

    # check whether gex was found
    assert np.array_equal(data.var["feature_types"].unique(), ["Gene Expression"]), "Output X should only contain Gene Expression vars."

    # check whether gene was found
    assert "PDE4DIP" in data.var_names, 'Output should contain gex column "PDE4DIP".'


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))