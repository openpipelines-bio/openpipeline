import sys
import pytest
from pathlib import Path
import mudata as md
import numpy as np
import scanpy as sc
from pandas.testing import assert_series_equal

## VIASH START
meta = {
    'executable': './target/docker/qc/calculate_qc_metrics/calculate_qc_metrics',
    'resources_dir': "./resources_test/pbmc_1k_protein_v3/",
    'config': './src/qc/calculate_qc_metrics/config.vsh.yaml',
    'cpus': 2
}
## VIASH END


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.fixture
def input_mudata(input_path):
    return md.read_h5mu(input_path)


@pytest.fixture
def mudata_w_na_boolean_colum(tmp_path, input_mudata):
    input_var = input_mudata.mod['rna'].var
    input_var["custom"] = np.nan

    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_mudata.write(new_input_path)
    return new_input_path

@pytest.fixture
def mudata_w_false_boolean_colum(tmp_path, input_mudata):
    input_var = input_mudata.mod['rna'].var
    input_var["custom"] = True

    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_mudata.write(new_input_path)
    return new_input_path


@pytest.fixture
def mudata_w_true_boolean_colum(tmp_path, input_mudata):
    input_var = input_mudata.mod['rna'].var
    input_var["custom"] = True

    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_mudata.write(new_input_path)
    return new_input_path

@pytest.fixture
def mudata_w_random_boolean_column(tmp_path, input_mudata):
    input_var = input_mudata.mod['rna'].var
    input_var["custom"] = np.random.choice([True, False], len(input_var), p=[0.8, 0.2])

    new_input_path = tmp_path / "input_with_custom_col.h5mu"
    input_mudata.write(new_input_path)
    return new_input_path

def test_add_qc(run_component, input_path):
    run_component([
        "--input", input_path,
        "--output", "foo.h5mu",
        "--modality", "rna",
        "--top_n_vars", "10,20,90",
        "--output_compression", "gzip"
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

@pytest.mark.parametrize("edited_input_mudata", ["mudata_w_true_boolean_colum",
                                                 "mudata_w_false_boolean_colum",
                                                 "mudata_w_na_boolean_colum",
                                                 "mudata_w_random_boolean_column"])
def test_calculcate_qc_var_qc_metrics(run_component, request, edited_input_mudata, tmp_path):
    input_path = request.getfixturevalue(edited_input_mudata)
    output_path = tmp_path / "foo.h5mu"

    args = [
        "--input", str(input_path),
        "--output", str(output_path),
        "--modality", "rna",
        "--top_n_vars", "10,20,90",
        "--var_qc_metrics", "custom",
    ]
    if edited_input_mudata == "mudata_w_na_boolean_colum":
        args.extend(["--var_qc_metrics_fill_na_value", "True"])

    run_component(args)
    assert output_path.is_file()
    data_with_qc = md.read(output_path)
    for qc_metric in ('pct_custom', 'total_counts_custom'):
        assert qc_metric in data_with_qc.mod['rna'].obs

def test_compare_scanpy(run_component,
                        mudata_w_random_boolean_column,
                        input_mudata,
                        tmp_path):
    
    output_path = tmp_path / "foo.h5mu"

    run_component([
        "--input", str(mudata_w_random_boolean_column),
        "--output", str(output_path),
        "--modality", "rna",
        "--top_n_vars", "10,20,90",
        "--var_qc_metrics", "custom",
        ])
    assert output_path.is_file()

    component_data = md.read(output_path)
    rna_mod = component_data.mod['rna']

    sc.pp.calculate_qc_metrics(
        input_mudata.mod['rna'],
        expr_type="counts",
        var_type="genes",
        qc_vars=["custom"],
        percent_top=[10,20,90],
        use_raw=False,
        inplace=True,
        log1p=False
    )
    scanpy_var =  input_mudata.mod['rna'].var
    component_var = rna_mod.var

    vars_to_compare = {
        'pct_dropout': 'pct_dropout_by_counts',
        'num_nonzero_obs': 'n_cells_by_counts',
        'obs_mean': 'mean_counts',
        'total_counts': 'total_counts'
    }
    for from_var, to_var in vars_to_compare.items():
        assert_series_equal(component_var[from_var],
                            scanpy_var[to_var],
                            check_names=False)


    scanpy_obs =  input_mudata.mod['rna'].obs
    component_obs = rna_mod.obs
    obs_to_compare = {
        'num_nonzero_vars': 'n_genes_by_counts',
        'pct_custom': 'pct_counts_custom',
        'total_counts_custom': 'total_counts_custom',
        'total_counts': 'total_counts'
    }
    obs_to_compare |= {f'pct_of_counts_in_top_{i}_vars': f'pct_counts_in_top_{i}_genes' 
                       for i in (10, 20, 90)}
    for from_obs, to_obs in obs_to_compare.items():
        assert_series_equal(component_obs[from_obs],
                            scanpy_obs[to_obs],
                            check_names=False)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
