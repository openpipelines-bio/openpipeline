import pytest
import pandas as pd
import uuid
from anndata import AnnData
from mudata import MuData, read_h5mu


@pytest.fixture
def h5mu():
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
        var = pd.DataFrame([["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
                           index=df.columns, columns=["Feat", "sample_id_var"])
        ad1 = AnnData(df, var=var)
        var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
        ad2 = AnnData(df, var=var2)
        return MuData({'mod1': ad1, 'mod2': ad2})


@pytest.fixture
def write_temp_h5mu(tmp_path):
    def wrapper(test_h5mu):
        test_h5mu_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
        test_h5mu.write_h5mu(test_h5mu_path)
        return test_h5mu_path
    return wrapper


def test_rename_var(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "output.h5mu"

    run_component(["--input", write_temp_h5mu(h5mu),
                   "--modality", "mod1",
                   "--var_key_input", "Feat",
                   "--var_key_output", "Foo",
                   "--output", output
                   ])

    assert output.is_file(), "Some output file must have been created."
    output_data = read_h5mu(output)

    assert 'Feat' not in output_data.mod['mod1'].var.columns
    assert 'Foo' in output_data.mod['mod1'].var.columns
    assert 'Foo' not in output_data.mod['mod2'].var.columns
    assert 'Foo' not in output_data.var.columns


def test_rename_var_without_key(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "output.h5mu"
    run_component(["--input", write_temp_h5mu(h5mu),
                   "--modality", "mod1",
                   "--var_key_output", "Foo",
                   "--output", output
                   ])

    assert output.is_file(), "Some output file must have been created."
    output_data = read_h5mu(output)

    assert 'Foo' in output_data.mod['mod1'].var.columns
    assert all(output_data.mod['mod1'].var.index == output_data.mod['mod1'].var["Foo"])


if __name__ == '__main__':
    exit(pytest.main([__file__]))
