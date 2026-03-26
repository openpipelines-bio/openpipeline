import sys
import os
import scanpy as sc
import pytest
import mudata as mu
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
reference_file = f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"


def log_normalize(adata):
    adata_lognorm = adata.copy()
    sc.pp.normalize_total(adata_lognorm, target_sum=1e4)
    sc.pp.log1p(adata_lognorm)
    adata.layers["log_normalized"] = adata_lognorm.X
    return adata


@pytest.fixture
def input_path(write_mudata_to_file):
    mdata = mu.read_h5mu(input_file)
    adata = mdata.mod["rna"].copy()
    adata_lognorm = log_normalize(adata)
    mdata.mod["rna"] = adata_lognorm
    return write_mudata_to_file(mdata)


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--input_var_gene_names",
            "gene_symbol",
            "--reference",
            reference_file,
            "--reference_var_input",
            "highly_variable",
            "--reference_obs_target",
            "cell_ontology_class",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == [
        "singler_pred",
        "singler_probability",
        "singler_delta_next",
        "singler_pruned_labels",
    ]
    assert list(output_mudata.mod["rna"].obsm.keys()) == ["singler_scores"]

    predictions = output_mudata.mod["rna"].obs["singler_pred"]
    assert predictions.dtype == "category", (
        "Calculated predictions should be category dtype"
    )
    assert predictions.notna().all(), "None of the predictions should be NA"

    pruned_predictions = output_mudata.mod["rna"].obs["singler_pruned_labels"]
    assert pruned_predictions.dtype == "category", (
        "Pruned predictions should be category dtype"
    )
    assert pruned_predictions.isna().any(), (
        "Some of the pruned predictions should be NA"
    )
    assert not all(pruned_predictions.isna()), "Not all pruned predictions should be NA"

    delta_next = output_mudata.mod["rna"].obs["singler_delta_next"]
    assert delta_next.dtype == "float", "Calculated delta next should be float dtype"

    scores = output_mudata.mod["rna"].obsm["singler_scores"]
    assert scores.dtype == "float", "Calculated scores should be float dtype"
    assert all(-1 <= value <= 1 for value in scores.flatten()), (
        ".obsm at singler_scores has values outside the range [-1, 1]"
    )


def test_params(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--reference",
            reference_file,
            "--reference_obs_target",
            "cell_ontology_class",
            "--input_reference_gene_overlap",
            "1000",
            "--reference_var_gene_names",
            "ensemblid",
            "--reference_var_input",
            "highly_variable",
            "de_n_genes",
            "600",
            "--de_method",
            "wilcox",
            "--quantile",
            "0.75",
            "--fine_tuning_threshold",
            "0.1",
            "--prune",
            "False",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == [
        "singler_pred",
        "singler_probability",
        "singler_delta_next",
    ]
    assert list(output_mudata.mod["rna"].obsm.keys()) == ["singler_scores"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
