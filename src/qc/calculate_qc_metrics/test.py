import sys
import pytest
from pathlib import Path
import mudata as md
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/qc/calculate_qc_metrics/calculate_qc_metrics',
    'resources_dir': "./resources_test/pbmc_1k_protein_v3/",
    'config': './src/qc/calculate_qc_metrics/config.vsh.yaml',
    'cpus': 2
}

@pytest.fixture
def viash_executable():
    return "./bin/viash"


## VIASH END
resources_dir = meta["resources_dir"]
input_sample_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_add_qc(run_component):
    run_component([
        "--input", input_sample_file,
        "--output", "foo.h5mu",
        "--modality", "rna",
        "--top_n_vars", "10,20,90"
        ])
    
    assert Path("foo.h5mu").is_file()
    data_with_qc = md.read("foo.h5mu")
    var, obs = data_with_qc.mod['rna'].var, data_with_qc.mod['rna'].obs
    for top_n_vars in ("10", "20", "90"):
        assert f"pct_of_counts_in_top_{top_n_vars}_vars" in obs
    assert "total_counts" in obs
    assert "num_nonzero_vars" in obs
    assert "pct_dropout" in var
    assert "num_nonzero_obs" in var
    assert "obs_mean" in var
    assert "total_counts" in var

def test_calculcate_qc_var_qc_metrics(run_component, tmp_path):
    input_data = md.read_h5mu(input_sample_file)
    input_var = input_data.mod['rna'].var
    input_var["custom"] = np.random.choice([True, False], len(input_var), p=[0.8, 0.2])
    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_data.write(new_input_path)
    run_component([
        "--input", new_input_path,
        "--output", "foo.h5mu",
        "--modality", "rna",
        "--top_n_vars", "10,20,90",
        "--var_qc_metrics", "custom",
        ])
    assert Path("foo.h5mu").is_file()
    data_with_qc = md.read("foo.h5mu")
    for qc_metric in ('pct_custom', 'total_counts_custom'):
        assert qc_metric in data_with_qc.mod['rna'].obs

if __name__ == "__main__":
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))
