import sys
import pytest
import mudata as mu
import anndata as ad
import numpy as np

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


@pytest.mark.parametrize(
    "n_obs,metric",
    [
        # scanpy uses pynndescent instead of sklearn's brute force KNeighborsTransformer
        # for >= 4096 observations with a non-euclidean metric, or >= 8192 observations.
        (5000, "cosine"),
        (9000, "euclidean"),
    ],
)
def test_run_pynndescent(run_component, tmp_path, n_obs, metric):
    rng = np.random.default_rng(0)
    rna = ad.AnnData(
        X=rng.poisson(1, size=(n_obs, 20)).astype(np.float32),
        obsm={"X_pca": rng.normal(size=(n_obs, 10)).astype(np.float32)},
    )
    input_path = tmp_path / "input.h5mu"
    mu.MuData({"rna": rna}).write_h5mu(input_path)
    output = tmp_path / "output.h5mu"
    num_neighbors = 15
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output),
            "--metric",
            metric,
            "--num_neighbors",
            str(num_neighbors),
        ]
    )

    assert output.is_file(), "No output was created."
    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "neighbors" in rna_out.uns
    for key in ["distances", "connectivities"]:
        assert key in rna_out.obsp, f"Output should have .obsp['{key}']"
        assert rna_out.obsp[key].shape == (n_obs, n_obs)
    # Every observation should have at least num_neighbors - 1 neighbors (itself excluded).
    # The pynndescent graph can also store explicit zeros (e.g. the observation itself),
    # so the number of stored entries per row is not exactly num_neighbors - 1.
    distances = rna_out.obsp["distances"].copy()
    distances.eliminate_zeros()
    assert (np.diff(distances.indptr) >= num_neighbors - 1).all()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
