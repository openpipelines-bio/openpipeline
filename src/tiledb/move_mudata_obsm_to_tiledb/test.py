import sys
import pytest
import boto3
import os
from moto.server import ThreadedMotoServer
import socket
import tiledbsoma
import subprocess
import re
import pandas as pd
import numpy as np
import mudata
import json


## VIASH START
meta = {
    "executable": "target/executable/tiledb/move_mudata_obsm_to_tiledb/move_mudata_obsm_to_tiledb",
    "resources_dir": "./resources_test",
    "cpus": 2,
    "config": "./src/tiledb/move_mudata_obsm_to_tiledb/config.vsh.yaml",
}
sys.path.append("src/utils")
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

input_dir = f"{meta['resources_dir']}/tiledb/pbmc_1k_protein_v3_mms"


@pytest.fixture
def input_mudata():
    return mudata.read_h5mu(
        f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
    )


@pytest.fixture
def input_mudata_extra_output_slot(input_mudata):
    new_obsm_key = pd.DataFrame(
        np.random.rand(input_mudata["rna"].n_obs, 5),
        index=input_mudata["rna"].obs_names,
        columns=pd.Index(["a", "b", "c", "d", "e"]),
    )
    input_mudata["rna"].obsm["test_input_slot"] = new_obsm_key
    return input_mudata


@pytest.fixture
def input_mudata_path(random_h5mu_path, input_mudata):
    output_path = random_h5mu_path()
    input_mudata.write(output_path)
    return output_path


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
    return moto_server


def test_key_already_exists_raises(
    run_component, input_mudata_path, initiated_database
):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_uri",
                "s3://test",
                "--endpoint",
                initiated_database,
                "--s3_region",
                "us-east-1",
                "--output_modality",
                "rna",
                "--input_mudata",
                str(input_mudata_path),
                "--modality",
                "rna",
                "--obsm_input",
                "X_leiden_harmony_umap",
            ]
        )
    assert re.search(
        r"ValueError: The following keys already exist in the database: X_leiden_harmony_umap",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_obsm_key_raises(run_component, initiated_database, input_mudata_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_uri",
                "s3://test",
                "--endpoint",
                initiated_database,
                "--s3_region",
                "us-east-1",
                "--output_modality",
                "rna",
                "--input_mudata",
                str(input_mudata_path),
                "--modality",
                "rna",
                "--obsm_input",
                "doesnotexist",
            ]
        )
    assert re.search(
        r"Not all \.obsm keys were found in the input!",
        err.value.stdout.decode("utf-8"),
    )


def test_add(
    run_component, initiated_database, input_mudata_extra_output_slot, random_h5mu_path
):
    input_path = random_h5mu_path()
    input_mudata_extra_output_slot.write(input_path)
    run_component(
        [
            "--input_uri",
            "s3://test",
            "--endpoint",
            initiated_database,
            "--s3_region",
            "us-east-1",
            "--output_modality",
            "rna",
            "--input_mudata",
            str(input_path),
            "--modality",
            "rna",
            "--obsm_input",
            "test_input_slot",
        ]
    )
    obsm_key_uri = "s3://test/ms/rna/obsm/test_input_slot"
    tiledb_config = {
        "vfs.s3.no_sign_request": "false",
        "vfs.s3.region": "us-east-1",
        "vfs.s3.endpoint_override": initiated_database,
    }
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)
    with tiledbsoma.open(uri=obsm_key_uri, mode="r", context=context) as open_array:
        obsm_data = open_array.read().coos().concat().to_scipy().todense()
        assert obsm_data.shape == (713, 5)
        original_data = (
            input_mudata_extra_output_slot["rna"].obsm["test_input_slot"].to_numpy()
        )
        np.testing.assert_allclose(original_data, obsm_data)
        assert json.loads(open_array.metadata["column_index"]) == [
            "a",
            "b",
            "c",
            "d",
            "e",
        ]


if __name__ == "__main__":
    sys.exit(pytest.main(["-s", __file__]))
