import sys
import pytest
import mudata as mu

## VIASH START
meta = {
    "executable": "./target/docker/graph/neighbors/find_neighbors",
    "name": "find_neighbors",
    "resources_dir": "resources_test/",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5mu"


def test_run(run_component, tmp_path):
    output = tmp_path / "output.h5mu"

    cmd_pars = [
        meta["executable"],
        "--input",
        input,
        "--output",
        str(output),
        "--obsm_input",
        "X_pca",
        "--uns_output",
        "foo_neigh",
        "--obsp_distances",
        "bar_dist",
        "--obsp_connectivities",
        "baz_conn",
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert "foo_neigh" in rna_out.uns, "Output should have .uns['foo_neigh']"
    assert "baz_conn" in rna_out.obsp, "Output should have .obsp['baz_conn']"
    assert "bar_dist" in rna_out.obsp, "Output should have .obsp['bar_dist']"
    assert "foo_neigh" not in rna_in.uns, "Output should not have .uns['foo_neigh']"
    assert "baz_conn" not in rna_in.obsp, "Input should not have .obsp['baz_conn']"
    assert "bar_dist" not in rna_in.obsp, "Input should not have .obsp['bar_dist']"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
