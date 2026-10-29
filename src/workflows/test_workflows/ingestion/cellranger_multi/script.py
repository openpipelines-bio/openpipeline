from mudata import read_h5mu
import sys
import pytest

##VIASH START
par = {"input": "input.h5mu", "input_filtered": None}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    for input_path in par["input"]:
        input_mudata = read_h5mu(input_path)

        assert list(input_mudata.mod.keys()) == ["rna", "prot", "vdj_t"]
        assert list(input_mudata.uns.keys()) == ["metrics_cellranger"]
        expected_metrics = [
            "Category",
            "Library Type",
            "Grouped By",
            "Group Name",
            "Metric Name",
            "Metric Value",
        ]
        assert (
            input_mudata.uns["metrics_cellranger"].columns.to_list() == expected_metrics
        )


def test_filtered_data_has_fewer_cells():
    if not par.get("input_filtered"):
        pytest.skip("No filtered input provided.")
    for raw_path, filtered_path in zip(par["input"], par["input_filtered"]):
        raw_mudata = read_h5mu(raw_path)
        filtered_mudata = read_h5mu(filtered_path)

        assert filtered_mudata.mod["rna"].n_obs < raw_mudata.mod["rna"].n_obs, (
            "Expected the filtered count matrix to contain fewer cells than the "
            "raw count matrix when --output_filtered_data is enabled."
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
