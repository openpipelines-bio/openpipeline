import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    'executable': './target/executable/filter/remove_modality/remove_modality',
    'resources_dir': './resources_test/',
    'cpus': 2
}
## VIASH END

input_sample_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_remove_component(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_sample_file,
        "--modality", "rna",
        "--output", str(output_path),
        "--output_compression", "gzip"
    ])
    assert output_path.is_file()
    output = read_h5mu(output_path)
    assert list(output.mod.keys()) == ["prot"]

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
