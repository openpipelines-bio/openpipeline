import sys
import pytest
import boto3
import os
from moto.server import ThreadedMotoServer
from contextlib import contextmanager
import socket
import tiledbsoma
import subprocess
import re
import pandas as pd
import numpy as np
import mudata
import requests


## VIASH START
meta = {
    "executable": "target/executable/tiledb/move_mudata_obs_to_tiledb/move_mudata_obs_to_tiledb",
    "resources_dir": "./resources_test",
    "cpus": 2,
    "config": "./src/tiledb/move_mudata_obs_to_tiledb/config.vsh.yaml",
}
sys.path.append("src/utils")
## VIASH END

sys.path.append(meta["resources_dir"])

input_dir = f"{meta['resources_dir']}/tiledb/pbmc_1k_protein_v3_mms"


@pytest.fixture
def input_mudata():
    return mudata.read_h5mu(
        f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
    )


@pytest.fixture
def input_mudata_extra_output_slot(input_mudata):
    new_obs_col = pd.Series(
        np.random.rand(input_mudata["rna"].n_obs, 1).ravel(),
        index=input_mudata["rna"].obs_names,
    )
    input_mudata["rna"].obs["test_input_slot"] = new_obs_col
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


@contextmanager
def managed_moto_server(*args, **kwargs):
    server = ThreadedMotoServer(*args, **kwargs)
    server.start()
    try:
        yield server
    finally:
        server.stop()


@pytest.fixture(scope="function")
def moto_server(aws_credentials):
    """Fixture to run a mocked AWS server for testing."""
    # Note: pass `port=0` to get a random free port.
    with managed_moto_server(
        ip_address=socket.gethostbyname(socket.gethostname()), port=0
    ) as moto_server:
        yield moto_server


@pytest.fixture
def initiated_database(moto_server):
    host, port = moto_server.get_host_and_port()
    server_uri = f"http://{host}:{port}"
    client = boto3.client("s3", endpoint_url=server_uri, region_name="us-east-1")
    client.create_bucket(Bucket="test")

    def raise_(ex):
        raise ex

    for root, _, files in os.walk(input_dir, onerror=raise_):
        for filename in files:
            local_path = os.path.join(root, filename)
            relative_path = os.path.relpath(local_path, input_dir)
            client.upload_file(local_path, "test", relative_path)
    client.close()
    yield server_uri
    requests.post(f"{server_uri}/moto-api/reset")


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
                "--obs_input",
                "filter_with_counts",
            ]
        )
    assert re.search(
        r"ValueError: The following keys already exist in the database: filter_with_counts",
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
                "--obs_input",
                "doesnotexist",
            ]
        )
    assert re.search(
        r"Not all \.obs keys were found in the input!",
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
            "--obs_input",
            "test_input_slot",
        ]
    )
    obs_key_uri = "s3://test/obs"
    tiledb_config = {
        "vfs.s3.no_sign_request": "false",
        "vfs.s3.region": "us-east-1",
        "vfs.s3.endpoint_override": initiated_database,
    }
    context = tiledbsoma.SOMATileDBContext(tiledb_config=tiledb_config)
    with tiledbsoma.open(uri=obs_key_uri, mode="r", context=context) as open_array:
        obs_data = open_array.read().concat().to_pandas()
        obs_data = obs_data.set_index("cell_id")
        assert "test_input_slot" in obs_data.columns
        mudata_col = mudata.read_h5ad(input_path, mod="rna").obs["test_input_slot"]
        pd.testing.assert_series_equal(
            obs_data.loc[:, "test_input_slot"],
            mudata_col,
            check_like=True,
            check_names=False,
        )


def test_output_folder(
    run_component,
    initiated_database,
    input_mudata_extra_output_slot,
    random_h5mu_path,
    tmp_path,
):
    input_path = random_h5mu_path()
    input_mudata_extra_output_slot.write(input_path)
    output_path = tmp_path / "tiledb_out"
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
            "--obs_input",
            "test_input_slot",
            "--output_tiledb",
            output_path,
        ]
    )
    assert output_path.is_dir()
    obs_key_uri = output_path / "obs"
    assert obs_key_uri.is_dir()
    print(list(obs_key_uri.iterdir()), file=sys.stderr, flush=True)
    with tiledbsoma.open(
        uri=f"file://{str(obs_key_uri.resolve())}",
        mode="r",
        context=tiledbsoma.SOMATileDBContext(),
    ) as open_array:
        obs_data = open_array.read().concat().to_pandas()
        obs_data = obs_data.set_index("cell_id")
        assert "test_input_slot" in obs_data.columns
        mudata_col = mudata.read_h5ad(input_path, mod="rna").obs["test_input_slot"]
        pd.testing.assert_series_equal(
            obs_data.loc[:, "test_input_slot"],
            mudata_col,
            check_like=True,
            check_names=False,
        )


if __name__ == "__main__":
    sys.exit(pytest.main(["-s", __file__]))
