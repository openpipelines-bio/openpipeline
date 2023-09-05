import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu

## VIASH START
meta = {
    'executable': './target/docker/metadata/join_uns_to_obs/join_uns_to_obs',
    'config': './src/metadata/join_uns_to_obs/config.vsh.yml'
}
## VIASH END

@pytest.fixture
def test_mudata(tmp_path):
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A", "sample1"], ["B", "sample2"]], index=df.index, columns=["Obs", "sample_id"])
    var = pd.DataFrame([["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
                        index=df.columns, columns=["Feat", "sample_id_var"])
    obsm = pd.DataFrame([["X", "W"]], index=pd.Index([0]), columns=["uns_col1", "uns_col2"])
    ad1 = AnnData(df, obs=obs, var=var, uns={"obsm1": obsm})
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)

    test_h5mu = tmp_path / "input.h5mu"
    mudata = MuData({'mod1': ad1, 'mod2': ad2})
    mudata.write_h5mu(test_h5mu)
    return test_h5mu

def test_join_uns_to_obs(test_mudata, run_component, tmp_path):
    output_file = tmp_path / "output.h5mu"
    run_component([
        "--input", str(test_mudata),
        "--modality", "mod1",
        "--uns_key", "obsm1",
        "--output", str(output_file)
    ])
    assert output_file.is_file()
    output_mudata = read_h5mu(output_file)
    assert 'obsm1' in output_mudata.mod['mod1'].uns
    pd.testing.assert_index_equal(output_mudata.mod['mod1'].obs.columns, 
                                  pd.Index(['Obs', 'sample_id', 'uns_col1', 'uns_col2']))
    pd.testing.assert_series_equal(output_mudata.mod['mod1'].obs['uns_col1'], 
                                   pd.Series(['X', 'X'], index=pd.Index(['obs1', 'obs2']), 
                                             name='uns_col1', dtype='category'))
    pd.testing.assert_series_equal(output_mudata.mod['mod1'].obs['uns_col2'], 
                                   pd.Series(['W', 'W'], index=pd.Index(['obs1', 'obs2']), 
                                             name='uns_col2', dtype='category'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
