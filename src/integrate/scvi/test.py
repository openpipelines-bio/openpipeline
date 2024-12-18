import pytest
from pathlib import Path
import scvi
import mudata as mu
from anndata.tests.helpers import assert_equal

## VIASH START
meta = {
    "executable": "./target/executable/integrate/scvi/scvi",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

import sys

sys.path.append(meta["resources_dir"])

from subset_vars import subset_vars
from set_var_index import set_var_index

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"
scvi_model_file = f"{meta['resources_dir']}/scvi"


@pytest.fixture
def create_scvi_model(random_path, tmp_path):
    def wrapper(reference_file):
        input_data = mu.read_h5mu(reference_file)
        adata = input_data.mod["rna"]
        input_modality = adata.copy()
        input_modality = set_var_index(input_modality, "ensemblid")

        input_modality = subset_vars(input_modality, subset_col="highly_variable")

        scvi.model.SCVI.setup_anndata(
            input_modality, batch_key="donor_id", labels_key="cell_ontology_class"
        )

        scvi_model = scvi.model.SCVI(
            input_modality,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=1,
        )
        scvi_model.train(max_epochs=10)

        # reference_data.write_h5mu(reference_file)
        scvi_model.save(scvi_model_file, save_anndata=True, overwrite=True)

        return scvi_model_file

    return wrapper


@pytest.fixture
def mudata_with_mod_rna_obs_batch(tmp_path, request):
    obs_batch, var_input, obsm_output, model = request.param

    if model:
        ref_model = create_scvi_model(reference_file)
    else:
        ref_model = None

    new_input_file = tmp_path / "input.h5mu"

    input_data = mu.read_h5mu(input_file)
    input_rna = input_data.mod["rna"]
    input_rna.obs[obs_batch] = "A"
    column_index = input_rna.obs.columns.get_indexer([obs_batch])
    input_rna.obs.iloc[slice(input_rna.n_obs // 2, None), column_index] = "B"
    input_data.write(new_input_file.name)

    return new_input_file.name, input_rna, obs_batch, var_input, obsm_output, ref_model


@pytest.mark.parametrize(
    "mudata_with_mod_rna_obs_batch",
    [("batch", None, None, None), ("batch2", "filter_with_hvg", "X_int", None)],
    indirect=True,
)
def test_scvi(run_component, mudata_with_mod_rna_obs_batch):
    new_input_file, input_rna, obs_batch, var_input, obsm_output, ref_model = (
        mudata_with_mod_rna_obs_batch
    )

    args = [
        "--input",
        new_input_file,
        "--modality",
        "rna",
        "--obs_batch",
        obs_batch,
        "--output",
        "output.h5mu",
        "--output_model",
        "test/",
        "--max_epochs",
        "1",
        "--n_obs_min_count",
        "10",
        "--n_var_min_count",
        "10",
        "--output_compression",
        "gzip",
    ]

    if var_input is not None:
        args.extend(["--var_input", var_input])
    if obsm_output is not None:
        args.extend(["--obsm_output", obsm_output])
    if ref_model is not None:
        args.extend(["--scvi_reference_model", ref_model])
    run_component(args)

    # check files
    assert Path("output.h5mu").is_file(), "Output file does not exist"
    assert Path("test").is_dir()
    assert Path("test/model.pt").is_file()

    # check output h5mu
    output_data = mu.read_h5mu("output.h5mu")
    output_rna = output_data.mod["rna"]
    assert (
        output_rna.n_obs == input_rna.n_obs
    ), f"Number of observations changed\noutput_data: {output_data}"
    assert (
        output_rna.n_vars == input_rna.n_vars
    ), f"Number of variables changed\noutput_data: {output_data}"
    assert len(output_rna.varm) == len(
        input_rna.varm
    ), f"Number of varm changed\noutput_data: {output_data}"

    expected_obsm_output = "X_scvi_integrated" if obsm_output is None else obsm_output
    assert (
        expected_obsm_output in output_rna.obsm
    ), f".obsm['{expected_obsm_output}'] not added\noutput_data: {output_data}"

    # assert that nothing else has changed
    del output_rna.obsm[expected_obsm_output]
    assert_equal(input_rna, output_rna)


def test_hvg_subsetting_helper():
    input_data = mu.read_h5mu(input_file)
    adata = input_data.mod["rna"]

    old_n_genes = adata.n_vars

    adata.var["highly_variable_features"] = False
    adata.var.iloc[
        : old_n_genes // 2, adata.var.columns.get_indexer(["highly_variable_features"])
    ] = True

    adata = subset_vars(adata, subset_col="highly_variable_features")

    # Correct number of genes is subsetted
    assert adata.n_vars == old_n_genes // 2
    # Only HVG are subsetted
    assert adata.var["highly_variable_features"].all()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
