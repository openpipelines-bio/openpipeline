import sys
import pytest
import boto3
import os
import socket
import tiledbsoma
import tiledbsoma.io
import mudata
from moto.server import ThreadedMotoServer


tiledbsoma.logging.debug()

## VIASH START
meta = {
    "executable": "target/executable/convert/from_tiledb_to_h5mu/from_tiledb_to_h5mu",
    "resources_dir": "./resources_test",
    "cpus": 2,
    "config": "./src/convert/from_tiledb_to_h5mu/config.vsh.yaml",
}
sys.path.append("src/utils")
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

input_dir = f"{meta['resources_dir']}/tiledb/pbmc_1k_protein_v3_mms"


@pytest.fixture(scope="module")
def aws_credentials():
    """Mocked AWS Credentials for moto."""
    os.environ["AWS_ACCESS_KEY_ID"] = "testing"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "testing"
    os.environ["AWS_DEFAULT_REGION"] = "us-east-1"


@pytest.fixture(scope="function")
def moto_server(aws_credentials):
    """Fixture to run a mocked AWS server for testing."""
    ip_addr = socket.gethostbyname(socket.gethostname())
    # Note: pass `port=0` to get a random free port.
    server = ThreadedMotoServer(ip_address=ip_addr, port=0)
    server.start()
    host, port = server.get_host_and_port()
    yield f"http://{host}:{port}"
    server.stop()


@pytest.fixture
def initiated_database(moto_server):
    client = boto3.client("s3", endpoint_url=moto_server, region_name="us-east-1")
    client.create_bucket(Bucket="test")

    def raise_(ex):
        raise ex

    for root, _, files in os.walk(input_dir, onerror=raise_):
        for filename in files:
            local_path = os.path.join(root, filename)
            relative_path = os.path.relpath(local_path, input_dir)
            client.upload_file(local_path, "test", relative_path)
    client.close()
    return moto_server


def test_convert(run_component, initiated_database, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input_uri",
            "s3://test",
            "--endpoint",
            initiated_database,
            "--s3_region",
            "us-east-1",
            "--input_modality",
            "rna",
            "--output_modality",
            "rna",
            "--input_layers",
            "log_normalized;raw",
            "--output",
            output,
        ]
    )
    output_mudata = mudata.read_h5mu(output)
    assert output_mudata.mod_names == ["rna"]
    assert output_mudata.shape == (713, 33538)
    rna_modality = output_mudata["rna"]
    assert set(rna_modality.layers.keys()) == {"raw", "log_normalized"}
    assert set(rna_modality.obsm.keys()) == {
        "X_leiden_harmony_umap",
        "X_pca",
        "X_pca_integrated",
        "X_umap",
        "knn_distances",
        "knn_indices",
    }
    assert set(rna_modality.varm.keys()) == {"pca_loadings"}
    assert set(rna_modality.var.columns.to_list()) == {
        "filter_with_hvg",
        "pct_dropout",
        "gene_symbol",
        "obs_mean",
        "num_nonzero_obs",
        "feature_types",
        "filter_with_counts",
        "total_counts",
        "genome",
    }
    assert set(rna_modality.obs.columns.to_list()) == {
        "scrublet_doublet_score",
        "pct_of_counts_in_top_50_vars",
        "pct_of_counts_in_top_200_vars",
        "sample_id",
        "harmony_integration_leiden_1.0",
        "filter_with_scrublet",
        "filter_with_counts",
        "total_counts",
        "pct_of_counts_in_top_100_vars",
        "num_nonzero_vars",
        "pct_of_counts_in_top_500_vars",
    }


if __name__ == "__main__":
    sys.exit(pytest.main(["-s", __file__]))
