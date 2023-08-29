from pathlib import Path
from mudata import read_h5mu
import sys
import pytest

## VIASH START
meta = {
    'executable': './target/docker/convert/from_cellranger_multi_to_h5mu/from_cellranger_multi_to_h5mu',
    'resources_dir': 'resources_test/',
    'config': 'src/convert/from_cellranger_multi_to_h5mu/config.vsh.yaml'

}
## VIASH END


resources_dir, executable = meta["resources_dir"], meta["executable"]
cellranger_multi_output = f"{resources_dir}/10x_5k_anticmv/processed/10x_5k_anticmv.cellranger_multi.output.output"

def test_cellranger_multi_basic(run_component):
    run_component(["--input", cellranger_multi_output,
                   "--output", "output.h5mu",
                   "--output_compression", "gzip"])
    assert Path("output.h5mu").is_file()
    converted_data = read_h5mu("output.h5mu")
    assert list(converted_data.mod.keys()) == ['rna', 'prot', 'vdj_t']
    assert list(converted_data.uns.keys()) == ['metrics_cellranger']
    assert converted_data.uns['metrics_cellranger'].columns.to_list() == ['Category', 'Library Type', 
                                                                          'Grouped By', 'Group Name', 
                                                                          'Metric Name', 'Metric Value']
    
def test_cellranger_multi_to_h5mu_crispr(run_component):
    run_component(["--input", f"{resources_dir}/10x_5k_lung_crispr/processed/10x_5k_lung_crispr.cellranger_multi.output.output",
                   "--output", "output2.h5mu",
                   "--output_compression", "gzip"])
    assert Path("output2.h5mu").is_file()
    converted_data = read_h5mu("output2.h5mu")
    assert list(converted_data.mod.keys()) == ['rna', 'gdo']
    assert list(converted_data.uns.keys()) == ['metrics_cellranger']
    assert 'perturbation_efficiencies_by_feature' in converted_data.mod['gdo'].uns
    assert 'perturbation_efficiencies_by_target' in converted_data.mod['gdo'].uns
    assert 'feature_reference' not in converted_data.mod['rna'].uns
    assert 'feature_reference' in converted_data.mod['gdo'].uns


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))