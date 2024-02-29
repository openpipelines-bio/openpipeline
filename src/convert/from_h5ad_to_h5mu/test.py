import sys
import pytest
import mudata as mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    'resources_dir': 'resources_test'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_run(run_component, tmp_path):
    mdata = mu.read_h5mu(input)

    tmp_rna = tmp_path / "rna.h5ad"
    tmp_prot = tmp_path / "prot.h5ad"
    mdata.mod["rna"].write_h5ad(tmp_rna)
    mdata.mod["prot"].write_h5ad(tmp_prot)

    tmp_output = tmp_path / "output.h5mu"

    cmd_pars = [
        "--modality", "rna",
        "--input", str(tmp_rna),
        "--modality", "prot",
        "--input", str(tmp_prot),
        "--output", str(tmp_output),
        "--output_compression", "gzip"
    ]
    run_component(cmd_pars)

    assert tmp_output.is_file(), "No output was created."

    mdata2 = mu.read_h5mu(tmp_output)

    assert list(mdata2.mod.keys()) == ["rna", "prot"]
    
    assert_annotation_objects_equal(mdata2.mod["rna"], str(tmp_rna))
    assert_annotation_objects_equal(mdata2.mod["prot"], str(tmp_prot))

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))