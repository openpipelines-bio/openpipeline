from pathlib import Path
import pytest
import sys

import mudata as md
import numpy as np
import scanpy as sc
import muon as mu

## VIASH START
meta = {
    'executable': './target/docker/transform/tfidf/tfidf',
    'resources_dir': "./resources_test/cellranger_atac_tiny_bcl/counts/",
    'config': './src/transform/tfidf/config.vsh.yaml',
    'cpus': 2
}
## VIASH END

@pytest.fixture
def synthetic_example():
    atac = sc.AnnData(np.array([
        [0, 0, 0],
        [1, 0, 1],
        [10, 0, 0],
        [100, 0, 1],
        [1000, 0, 0]
    ]))
    atac.obs_names = ["A", "B", "C", "D", "E"]
    atac.var_names = ["x", "y", "z"]

    return md.MuData({"atac": atac})

@pytest.fixture
def example_mudata(tmp_path, synthetic_example):
    mdata_path = tmp_path / "example.h5mu"
    synthetic_example.write(mdata_path)

    return mdata_path

@pytest.fixture
def example_mudata_with_layer(tmp_path, synthetic_example):
    synthetic_example.mod["atac"].layers["atac_counts"] = synthetic_example.mod["atac"].X.copy()
    synthetic_example.mod["atac"].X = np.random.normal(size=synthetic_example.mod["atac"].X.shape)
    mdata_path = tmp_path / "example.h5mu"
    synthetic_example.write(mdata_path)

    return mdata_path

@pytest.fixture
def neurips_mudata(tmp_path):
    """From the `NeurIPS Multimodal Single-Cell Integration Challenge
    <https://www.kaggle.com/competitions/open-problems-multimodal/data>`
    
    Link is taken from the Moscot repository: 
    https://github.com/theislab/moscot/blob/cb53435c80fafe58046ead3c42a767fd0b818aaa/src/moscot/datasets.py#L67
    """
    adata = sc.read("../data/neurips_data.h5ad", backup_url="https://figshare.com/ndownloader/files/37993503")

    mdata = md.MuData({"atac": adata})
    mdata_path = tmp_path / "neurips.h5mu"
    mdata.write(mdata_path)

    return mdata_path

@pytest.fixture
def tiny_atac_mudata(tmp_path):
    resources_dir = Path(meta["resources_dir"])

    mdata = mu.read_10x_h5(resources_dir / "counts" / "filtered_peak_bc_matrix.h5")
    mdata_path = tmp_path / "tiny_atac.h5mu"
    mdata.write(mdata_path)

    return mdata_path

@pytest.mark.parametrize("mudata", ["example_mudata", "neurips_mudata", "tiny_atac_mudata"])
def test_output_layer(run_component, request, mudata, tmp_path):
    input_path = request.getfixturevalue(mudata)
    output_path = tmp_path / "foo.h5mu"

    args = [
        "--input", str(input_path),
        "--output", str(output_path),
        "--modality", "atac",
    ]

    run_component(args)
    assert output_path.is_file()
    output_mdata = md.read(output_path)

    assert "tfidf" in output_mdata.mod["atac"].layers.keys()

@pytest.mark.parametrize("mudata", ["example_mudata"])
def test_calculations_correctness(request, run_component, mudata, tmp_path):
    input_path = request.getfixturevalue(mudata)
    output_path = tmp_path / "foo.h5mu"

    args = [
        "--input", str(input_path),
        "--output", str(output_path),
        "--modality", "atac",
    ]

    run_component(args + ["--scale_factor", "10000", "--output_layer", "tfidf_10000"])
    assert output_path.is_file()
    output_mdata = md.read(output_path)

    assert np.allclose(
        output_mdata.mod["atac"].layers["tfidf_10000"].toarray(),
        np.array([[    np.nan, np.nan,      np.nan],
                  [0.0382461 , 0.    , 10.67027475],
                  [0.04135813, 0.    , 0.         ],
                  [0.04131346, 0.    , 5.7693107  ],
                  [0.04135813, 0.    , 0.         ]]),
        equal_nan=True
        )

    run_component(args + ["--scale_factor", "100", "--output_layer", "tfidf_100"])
    output_mdata = md.read(output_path)

    assert np.allclose(
        output_mdata.mod["atac"].layers["tfidf_100"].toarray(),
        np.array([[    np.nan, np.nan,     np.nan],
                  [0.01765529, 0.    , 4.92564555],
                  [0.02072352, 0.    , 0.        ],
                  [0.02067929, 0.    , 0.86213192],
                  [0.02072352, 0.    , 0.        ]]),
        equal_nan=True
    )
    

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
