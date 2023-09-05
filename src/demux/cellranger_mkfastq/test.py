import sys
import pytest

## VIASH START
meta = {
    "functionality_name": "cellranger_mkfastq",
    "resources_dir": "resources_test"
}
## VIASH END

input = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl"
sample_sheet = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl/sample_sheet.csv"

def test_run(run_component, tmp_path):
    output = tmp_path / "output"

    cmd_pars = [
        "--input", input,
        "--sample_sheet", sample_sheet,
        "--output", str(output)
    ]
    if meta['cpus']:
        cmd_pars.extend(["---cpus", str(meta['cpus'])])
    if meta['memory_gb']:
        cmd_pars.extend(["---memory", f"{meta['memory_gb']}GB"])

    run_component(cmd_pars)

    expected_dir: Path = output / "H35KCBCXY" / "test_sample"
    assert expected_dir.is_dir(), "No output was created."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))