import sys
import pytest
from mudata import read_h5mu, MuData
from anndata import AnnData
import pandas as pd
import uuid

## VIASH START
meta = {
    'executable': './target/docker/filter/intersect_obs/intersect_obs',
    'resources_dir': './resources_test/',
    'cpus': 2,
    'config': './src/filter/intersect_modalities/config.vsh.yaml'
}
## VIASH END

input_sample_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def generate_h5mu():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]],
                        index=df.columns, columns=["Feat"])
    ad1 = AnnData(df, obs=obs, var=var)
    df2 = pd.DataFrame([[7, 8, 9], [10, 11, 12]], index=["obs2", "obs3"], columns=df.columns)
    var2 = pd.DataFrame(["d", "e", "g"], index=df2.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df2.index, columns=["Obs"])
    ad2 = AnnData(df2, obs=obs2, var=var2)
    tmp_mudata = MuData({'mod1': ad1, 'mod2': ad2})
    return tmp_mudata

@pytest.fixture
def sample_mudata(generate_h5mu, tmp_path):
    filename = f"{uuid.uuid4()}.h5mu"
    output_file = tmp_path / filename
    generate_h5mu.write(output_file)
    return output_file
    

def test_intersect_obs(run_component, sample_mudata, tmp_path):
    output_path = tmp_path / f"{uuid.uuid4()}.h5mu"

    # run component
    run_component([
        "--input", sample_mudata,
        "--modalities", "mod1:mod2",
        "--output", str(output_path),
        "--output_compression", "gzip"
    ])
    assert output_path.is_file()
    output = read_h5mu(output_path)
    assert list(output.mod.keys()) == ["mod1", "mod2"]
    assert output.obs_names.tolist() == ["obs2"]

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))