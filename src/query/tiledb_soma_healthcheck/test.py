import sys
import pytest

## VIASH START
meta = {
    "resources_dir": "./resources_test/",
    "executable": "./target/executable/query/tiledb_hea",
    "config": "/home/di/code/openpipeline/src/query/cellxgene_census/config.vsh.yaml",
}
## VIASH END


def test_healthcheck(run_component):
    run_component(
        [
            "--input_uri",
            "s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/",
            "--s3_region",
            "us-west-2",
        ]
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
