import subprocess
import sys

import h5py
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

## VIASH START
meta = {
    "executable": "./target/executable/process_10xh5/filter_10xh5/filter_10xh5",
    "resources_dir": "resources_test/pbmc_1k_protein_v3",
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"


def read_matrix(path):
    with h5py.File(path, "r") as h5:
        shape = tuple(h5["matrix/shape"][()])
        mat = sp.csc_matrix(
            (h5["matrix/data"][()], h5["matrix/indices"][()], h5["matrix/indptr"][()]),
            shape=shape,
        )
        barcodes = h5["matrix/barcodes"].asstr()[()]
        ids = h5["matrix/features/id"].asstr()[()]
        feature_types = h5["matrix/features/feature_type"].asstr()[()]
    return mat, barcodes, ids, feature_types


def run_component(output_path, *args):
    subprocess.run(
        [
            meta["executable"],
            "--input",
            input_path,
            "--output",
            str(output_path),
            "--verbose",
            *args,
        ],
        check=True,
    )


def test_filter(tmp_path):
    output_path = tmp_path / "out.h5"
    run_component(
        output_path, "--min_library_size", "1000", "--min_cells_per_gene", "300"
    )
    assert output_path.is_file()

    orig_mat, orig_barcodes, orig_ids, _ = read_matrix(input_path)
    mat, barcodes, ids, _ = read_matrix(output_path)

    # compute expected result independently
    cell_mask = np.asarray(orig_mat.sum(axis=0)).ravel() >= 1000
    sub = orig_mat[:, cell_mask]
    feature_mask = np.asarray((sub > 0).sum(axis=1)).ravel() >= 300
    expected = sub[feature_mask, :]

    np.testing.assert_array_equal(barcodes, orig_barcodes[cell_mask])
    np.testing.assert_array_equal(ids, orig_ids[feature_mask])
    assert mat.shape == expected.shape
    assert (mat != expected).nnz == 0

    # filtering was effective but did not remove everything
    assert 100 <= mat.shape[1] < orig_mat.shape[1]
    assert 100 <= mat.shape[0] < orig_mat.shape[0]

    # output must be readable by scanpy (as a 10x h5 file)
    adata = sc.read_10x_h5(output_path, gex_only=False)
    assert adata.shape == (mat.shape[1], mat.shape[0])

    # dtypes, file-level attributes and other feature metadata are kept
    with (
        h5py.File(input_path, "r") as h5_in,
        h5py.File(output_path, "r") as h5_out,
    ):
        assert set(h5_out.attrs.keys()) == set(h5_in.attrs.keys())
        for key in h5_in.attrs:
            np.testing.assert_array_equal(h5_out.attrs[key], h5_in.attrs[key])
        assert set(h5_out["matrix/features"].keys()) == set(
            h5_in["matrix/features"].keys()
        )
        np.testing.assert_array_equal(
            h5_out["matrix/features/_all_tag_keys"][()],
            h5_in["matrix/features/_all_tag_keys"][()],
        )
        for name in ["data", "indices", "indptr", "shape", "barcodes"]:
            assert h5_out[f"matrix/{name}"].dtype == h5_in[f"matrix/{name}"].dtype


def test_keep_feature_types(tmp_path):
    output_path = tmp_path / "out.h5"
    run_component(
        output_path,
        "--min_cells_per_gene",
        "100000",
        "--keep_feature_types",
        "Antibody Capture",
    )
    _, _, _, orig_feature_types = read_matrix(input_path)
    mat, barcodes, _, feature_types = read_matrix(output_path)

    n_adt = int((orig_feature_types == "Antibody Capture").sum())
    assert n_adt > 0
    assert mat.shape[0] == n_adt
    assert set(feature_types) == {"Antibody Capture"}


def test_no_filtering(tmp_path):
    output_path = tmp_path / "out.h5"
    run_component(output_path)
    orig_mat, orig_barcodes, orig_ids, _ = read_matrix(input_path)
    mat, barcodes, ids, _ = read_matrix(output_path)
    assert (mat != orig_mat).nnz == 0
    np.testing.assert_array_equal(barcodes, orig_barcodes)
    np.testing.assert_array_equal(ids, orig_ids)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
