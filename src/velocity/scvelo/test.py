import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    "name": "./target/executable/projection/scvelo/scvelo",
    "resources_dir": "./resources_test/",
}
## VIASH END

input_loom = f"{meta['resources_dir']}/cellranger_tiny.loom"


def test_scvelo(run_component, tmp_path):
    output_dir = tmp_path / "foo"
    run_component(
        [
            "--input",
            input_loom,
            "--output",
            str(output_dir),
            "--output_compression",
            "gzip",
        ]
    )

    assert output_dir.is_dir()
    assert (output_dir / "scvelo_proportions.pdf").is_file()
    assert (output_dir / "scvelo_embedding.pdf").is_file()
    assert (output_dir / "scvelo_graph.pdf").is_file()
    assert (output_dir / "proportions.txt").is_file()
    assert (output_dir / "foo.h5mu").is_file()

    output_data = read_h5mu(output_dir / "foo.h5mu")
    assert "rna_velocity" in output_data.mod.keys()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
