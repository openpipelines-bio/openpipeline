import sys
import pytest
import mudata as mu
from pathlib import Path
import pandas as pd

## VIASH START
meta = {
    "executable": "./target/executable/compression/compress_h5mu/compress_h5mu",
    "resources_dir": "resources_test/concat_test_data/",
    "config": "src/compression/compress_h5mu/config.vsh.yaml",
}
## VIASH END


input_file = Path(
    f"{meta['resources_dir']}/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"
)


def compare_anndata(first, second):
    for attr_name in ("obs", "var"):
        pd.testing.assert_frame_equal(
            getattr(first, attr_name), getattr(second, attr_name)
        )


@pytest.mark.parametrize("compression_type", ["gzip", "lzf"])
def test_compress_h5mu(run_component, tmp_path, compression_type):
    output_file = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_file),
            "--output",
            str(output_file),
            "--compression",
            compression_type,
        ]
    )

    # check whether file exists
    assert output_file.is_file(), "Output file does not exist"

    # read output mudata
    output = mu.read_h5mu(output_file)
    uncompressed_h5mu = mu.read_h5mu(input_file)
    for attr_name in ("obs", "var"):
        pd.testing.assert_frame_equal(
            getattr(output, attr_name), getattr(uncompressed_h5mu, attr_name)
        )
    for mod_name in uncompressed_h5mu.mod:
        assert (
            mod_name in output.mod
        ), f"{mod_name} found in uncompressed file, but not in compressed output file."
        mod_compressed = output.mod[mod_name]
        mod_uncompressed = uncompressed_h5mu.mod[mod_name]
        compare_anndata(mod_compressed, mod_uncompressed)
    assert output_file.stat().st_size < input_file.stat().st_size


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
