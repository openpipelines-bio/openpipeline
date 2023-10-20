


import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
from pathlib import Path
from subprocess import CalledProcessError

## VIASH START
meta = {
    'executable': './target/docker/metadata/grep_annotation_column/grep_annotation_column',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2,
    'config': '/home/di/code/openpipelines-multisample/src/metadata/grep_annotation_column/config.vsh.yaml'
}
## VIASH END

@pytest.fixture
def generate_h5mu(tmp_path):
    tmp_file = tmp_path / "test.h5mu"

    # generate data
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]],
                        index=df.columns, columns=["Feat"])
    ad1 = AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)
    tmp_mudata = MuData({'mod1': ad1, 'mod2': ad2})
    
    # write to file
    tmp_mudata.write_h5mu(tmp_file)
    return tmp_file

@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column(run_component, generate_h5mu, tmp_path, compression_format):
    output_path = tmp_path / "with_id.h5mu"
    
    # run component
    run_component([
        "--input", str(generate_h5mu),
        "--output", str(output_path),
        "--modality", "mod1",
        "--matrix", "var",
        "--input_column", "Feat",
        "--regex_pattern", "^var1",
        "--output_match_column", "test",
        "--output_compression", compression_format
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod['mod1'].var.columns.to_list()
    assert output_data.mod['mod1'].var['test'].to_list() == [True, False, False]

@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column_default(run_component, generate_h5mu, tmp_path, compression_format):
    output_path = tmp_path / "with_id.h5mu"
    
    # run component
    run_component([
        "--input", str(generate_h5mu),
        "--output", str(output_path),
        "--modality", "mod1",
        "--matrix", "var",
        "--regex_pattern", "^var1",
        "--output_match_column", "test",
        "--output_compression", compression_format
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod['mod1'].var.columns.to_list()
    assert output_data.mod['mod1'].var['test'].to_list() == [True, False, False]

@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column(run_component, generate_h5mu, tmp_path, compression_format):
    output_path = tmp_path / "with_id.h5mu"
    
    # run component
    run_component([
        "--input", str(generate_h5mu),
        "--output", str(output_path),
        "--modality", "mod1",
        "--matrix", "var",
        "--input_column", "Feat",
        "--regex_pattern", "^a",
        "--output_match_column", "test",
        "--output_compression", compression_format
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod['mod1'].var.columns.to_list()
    assert output_data.mod['mod1'].var['test'].to_list() == [True, False, False]


@pytest.mark.parametrize("compression_format", ["gzip", "lzf"])
def test_grep_column_fraction_column(run_component, generate_h5mu, tmp_path, compression_format):
    output_path = tmp_path / "with_id.h5mu"
    
    # run component
    run_component([
        "--input", str(generate_h5mu),
        "--output", str(output_path),
        "--modality", "mod1",
        "--matrix", "var",
        "--input_column", "Feat",
        "--regex_pattern", "^a",
        "--output_match_column", "test",
        "--output_fraction_column", "test_output_fraction",
        "--output_compression", compression_format
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test" in output_data.mod['mod1'].var.columns.to_list()
    assert output_data.mod['mod1'].var['test'].to_list() == [True, False, False]
    output_data.mod['mod1'].obs['test_output_fraction'].to_list() == [1/6, 4/15]

def test_missing_column(run_component, generate_h5mu, tmp_path,):
    output_path = tmp_path / "with_id.h5mu"

    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", str(generate_h5mu),
            "--output", str(output_path),
            "--modality", "mod1",
            "--matrix", "var",
            "--input_column", "filliberke",
            "--regex_pattern", "^a",
            "--output_match_column", "test",
            "--output_compression", "gzip",
        ])
    assert "ValueError: Column filliberke could not be found for modality mod1" in \
        err.value.stdout.decode('utf-8')

def test_invalid_regex_pattern(run_component, generate_h5mu, tmp_path,):
    output_path = tmp_path / "with_id.h5mu"

    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", str(generate_h5mu),
            "--output", str(output_path),
            "--modality", "mod1",
            "--matrix", "var",
            "--input_column", "Feat",
            "--regex_pattern", "(a",
            "--output_match_column", "test",
            "--output_compression", "gzip",
        ])
    assert "ValueError: (a is not a valid regular expression pattern." in \
        err.value.stdout.decode('utf-8')


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-s"]))

