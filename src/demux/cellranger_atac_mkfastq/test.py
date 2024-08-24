import sys
from pathlib import Path
import pytest

## VIASH START
meta = {
    "functionality_name": "cellranger_mkfastq",
    "resources_dir": "resources_test"
}
## VIASH END

input_dir = Path(meta["resources_dir"]) / "cellranger_atac_tiny_bcl/bcl"
sample_sheet = Path(meta["resources_dir"]) / "cellranger_atac_tiny_bcl/bcl/layout.csv"

def test_run(run_component, tmp_path):
    output = tmp_path / "output"

    print("Input dir exists: ", input_dir.is_dir())
    print("Input dir content: ", list(input_dir.glob("*")))
    print("Sample sheet exists: ", sample_sheet.is_file())

    cmd_pars = [
        "--input", str(input_dir),
        "--csv", str(sample_sheet),
        "--output", str(output)
    ]

    run_component(cmd_pars)

    expected_dir: Path = output / "HJN3KBCX2" / "test_sample"
    assert expected_dir.is_dir(), "No output was created."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))