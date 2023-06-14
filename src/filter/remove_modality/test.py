import pytest
from pathlib import Path
from mudata import read_h5mu
from sys import exit

## VIASH START
meta = {
    'executable': './target/docker/filter/remove_modality/remove_modality',
    'resources_dir': './resources_test/',
    'cpus': 2
}
## VIASH END


meta['cpus'] = 1 if not meta['cpus'] else meta['cpus']

resources_dir = meta["resources_dir"]
input_sample_file = f"{resources_dir}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_remove_component(run_component):
    run_component(["--input", input_sample_file,
                   "--modality", "rna",
                   "--output", "removed_rna.h5mu",
                   "--output_compression", "gzip"])
    assert Path("removed_rna.h5mu").is_file()
    output = read_h5mu("removed_rna.h5mu")
    assert list(output.mod.keys()) == ["prot"]

if __name__ == '__main__':
    exit(pytest.main([__file__], plugins=["viashpy"]))
