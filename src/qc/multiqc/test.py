import sys
import pytest
import shutil

## VIASH START
## VIASH END

input_fastqc = meta["resources_dir"] + "/fastqc/"


def test_multiqc(run_component, tmp_path):
    output_path = tmp_path / "output"

    run_component(["--input", input_fastqc, "--output", str(output_path)])

    assert output_path.exists()
    assert (output_path / "multiqc_report.html").is_file()


def test_multiple_inputs(run_component, tmp_path):
    output_path = tmp_path / "output"
    dir1 = tmp_path / "dir1"
    dir2 = tmp_path / "dir2"

    # copy input to tmp_path
    shutil.copytree(input_fastqc, dir1)
    shutil.copytree(input_fastqc, dir2)

    run_component(
        ["--input", str(dir1), "--input", str(dir2), "--output", str(output_path)]
    )

    assert output_path.exists()
    assert (output_path / "multiqc_report.html").is_file()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
