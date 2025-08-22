import sys
import os
import pytest
import mudata as mu

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"


def test_simple_execution(run_component):
    output_file = "output.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--reference",
            reference_file,
            "--output",
            "output.h5mu",
            "--methods",
            "rf;svm",
        ]
    )

    # check whether file exists
    assert os.path.exists(output_file), "Output file does not exist"

    # read output mudata
    output = mu.read_h5mu(output_file)

    # check output
    expected_rna_obs_cols = ["popv_prediction"]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, (
            f"could not find columns .mod['rna'].obs['{col}']"
        )

    print(f"output: {output}", flush=True)


def test_popv_with_other_layer(run_component, tmp_path):
    input_h5mu = mu.read(input_file)
    input_h5mu.mod["rna"].layers["test"] = input_h5mu.mod["rna"].X.copy()
    input_h5mu.write_h5mu(tmp_path / "input.h5mu")
    run_component(
        [
            "--input",
            tmp_path / "input.h5mu",
            "--reference",
            reference_file,
            "--output",
            "output.h5mu",
            "--methods",
            "rf;svm;knn_on_scanorama;knn_on_scvi",
        ]
    )


def test_popv_with_non_overlapping_cells(run_component, tmp_path):
    input_h5mu = mu.read(input_file)

    # copy previous modalities
    rna_ad = input_h5mu.mod["rna"].copy()
    prot_ad = input_h5mu.mod["prot"].copy()

    # change obs_names such that the cells do not overlap
    rna_ad.obs_names = [f"rna_{x}" for x in rna_ad.obs_names]
    prot_ad.obs_names = [f"prot_{x}" for x in prot_ad.obs_names]

    # write new h5mu to file
    new_h5mu = mu.MuData({"rna": rna_ad, "prot": prot_ad})
    new_h5mu.write_h5mu(tmp_path / "input.h5mu")

    # run component
    run_component(
        [
            "--input",
            tmp_path / "input.h5mu",
            "--reference",
            reference_file,
            "--output",
            "output.h5mu",
            "--methods",
            "rf;svm;knn_on_scanorama",
        ]
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
