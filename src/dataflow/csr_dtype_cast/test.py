import sys
import pytest
from pathlib import Path
import mudata as mu
import anndata as ad
import numpy as np
import pandas as pd
import uuid
from scipy.sparse import random, csr_array


@pytest.fixture
def input_mdata():
    rng = np.random.default_rng(seed=1)
    random_counts = random(
        50000, 100,
        density=0.8,
        format='csr',
        dtype=np.uint32,
        random_state=rng
        )
    bad_dtype = random_counts.astype(np.float32)
    bad_dtype.indptr = bad_dtype.indptr.astype(np.uint32)
    bad_dtype.indices = bad_dtype.indices.astype(np.int64)
    del random_counts 
    mod1 = ad.AnnData(
        X=bad_dtype,
        obs=pd.DataFrame(index=pd.RangeIndex(50000)),
        var=pd.DataFrame(index=pd.RangeIndex(100))
        )

    return mu.MuData({"mod1": mod1})


@pytest.fixture
def input_mdata_path(tmp_path, input_mdata):
    output_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
    input_mdata.write(output_path)
    return output_path


def test_dtype_casting(
        run_component,
        input_mdata_path
        ):

    args = [
        "--input", input_mdata_path,
        "--output", "foo.h5mu",
        "--modality", "mod1",
        "--output_compression", "gzip",
        "--dtype", "int32"
    ]

    run_component(args)
    assert Path("foo.h5mu").is_file()
    mdata = mu.read("foo.h5mu")

    assert mdata.mod["mod1"].X.indices.dtype == np.int32, "Indices dtype not casted correctly"
    assert mdata.mod["mod1"].X.indptr.dtype == np.int32, "Indptr dtype not casted correctly"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
