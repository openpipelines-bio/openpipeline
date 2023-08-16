import pytest
import re
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
from subprocess import CalledProcessError

## VIASH START
meta = {
    'functionality_name': 'move_obsm_to_obs',
    'resources_dir': 'resources_test/',
    'executable': 'target/docker/metadata/move_obsm_to_obs/move_obsm_to_obs',
    'config': '/home/di/code/openpipeline/src/metadata/move_obsm_to_obs/config.vsh.yaml'
}
## VIASH END

@pytest.fixture
def h5mu():
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
        obs = pd.DataFrame([["A", "sample1"], ["B", "sample2"]], index=df.index, columns=["Obs", "sample_id"])
        var = pd.DataFrame([["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
                           index=df.columns, columns=["Feat", "sample_id_var"])
        obsm = {"obsm_key": pd.DataFrame([["foo", "bar"], ["lorem", "ipsum"]],
                                         index=obs.index, columns=["obsm_col1", "obsm_col2"])}
        ad1 = AnnData(df, obs=obs, var=var, obsm=obsm)
        var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
        obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
        ad2 = AnnData(df, obs=obs2, var=var2)
        return MuData({'mod1': ad1, 'mod2': ad2})

@pytest.fixture
def write_temp_h5mu(tmp_path):
    def wrapper(test_h5mu): 
        test_h5mu_path = tmp_path / "input.h5mu"
        test_h5mu.write_h5mu(test_h5mu_path.name)
        return test_h5mu_path.name
    return wrapper

@pytest.fixture
def h5mu_with_non_overlapping_observations(h5mu):
    h5mu.mod['mod1'].obsm['obsm_key'].index = pd.Index(["obs1", "doesnotexist"])
    return h5mu


def test_move_obsm_to_obs(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path/ "output.h5mu"
    run_component(["--input", write_temp_h5mu(h5mu),
                   "--modality", "mod1",
                   "--obsm_key", "obsm_key",
                   "--output", output
                   ])
    assert output.is_file(), "Some output file must have been created."
    output_data = read_h5mu(output)
    pd.testing.assert_index_equal(output_data.mod['mod1'].obs.index, pd.Index(['obs1', 'obs2']))
    pd.testing.assert_index_equal(output_data.mod['mod1'].obs.columns, 
                                  pd.Index(['Obs', 'sample_id', 'obsm_key_obsm_col1', 'obsm_key_obsm_col2']))
    assert 'obsm_key' not in output_data.mod['mod1'].obsm

def test_move_obsm_to_obs_non_overlapping_obs_fails(run_component, write_temp_h5mu,
                                                    h5mu_with_non_overlapping_observations, tmp_path):
    output = tmp_path/ "output.h5mu"
    # Mudata seems to handle this error, but keep this test in just in case mudata drops the ball.
    with pytest.raises(CalledProcessError) as err:
        run_component(["--input", write_temp_h5mu(h5mu_with_non_overlapping_observations),
                    "--modality", "mod1",
                    "--obsm_key", "obsm_key",
                    "--output", output
                    ])
    re.search(r"value.index does not match parent’s axis 0 names",
        err.value.stdout.decode('utf-8'))


def test_error_non_existing_modality(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path/ "output.h5mu"
    with pytest.raises(CalledProcessError) as err:
        run_component(["--input", write_temp_h5mu(h5mu),
                    "--modality", "foo",
                    "--obsm_key", "obsm_key",
                    "--output", output
                    ])
    re.search(r"ValueError: Modality foo does not exist\.",
        err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    exit(pytest.main([__file__]))