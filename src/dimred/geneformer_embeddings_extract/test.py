import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/dimred/geneformer_embeddings_extract/geneformer_embeddings_extract",
    "resources_dir": "resources_test",
    "config": "./src/dimred/geneformer_embeddings_extract/config.vsh.yaml",
}
## VIASH END

model_dir = f"{meta['resources_dir']}/geneformer/Geneformer-V1-10M"

# Geneformer-V1-10M: hidden size 256, input size 2048, no special tokens
EMB_DIMS = 256
PAD_TOKEN = 0


def make_tokens(n_cells, rng):
    """Token sequences of different lengths, so that the length sort inside
    EmbExtractor reorders the cells."""
    tokens = np.full((n_cells, 2048), PAD_TOKEN, dtype=np.uint16)
    lengths = rng.integers(20, 300, size=n_cells)
    for row, length in enumerate(lengths):
        # V1 dictionary: 0 = <pad>, 1 = <mask>, genes from 2 onwards
        tokens[row, :length] = rng.choice(
            np.arange(2, 20_000), size=length, replace=False
        )
    return tokens


@pytest.fixture
def tokenized_h5mu():
    rng = np.random.default_rng(0)
    n_cells = 24
    rna = ad.AnnData(
        X=np.zeros((n_cells, 3)),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    rna.obsm["geneformer_tokens"] = make_tokens(n_cells, rng)
    rna.uns["geneformer_tokenize"] = {
        "model_version": "V1",
        "model_input_size": 2048,
        "special_token": False,
        "pad_token": PAD_TOKEN,
        "cls_token": -1,
        "eos_token": -1,
    }
    prot = ad.AnnData(X=np.ones((n_cells, 2)), obs=pd.DataFrame(index=rna.obs_names))
    return mu.MuData({"rna": rna, "prot": prot})


@pytest.fixture
def input_path(tokenized_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    tokenized_h5mu.write_h5mu(path)
    return path


def run_embed(run_component, input_path, output_path, *extra):
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--model",
            model_dir,
            "--model_version",
            "V1",
            "--forward_batch_size",
            "5",
            *extra,
        ]
    )
    return mu.read_h5mu(output_path)


def test_embeddings_keep_input_order(
    run_component, input_path, tokenized_h5mu, random_h5mu_path
):
    output = run_embed(run_component, input_path, random_h5mu_path())
    embeddings = output.mod["rna"].obsm["X_geneformer"]
    assert embeddings.shape == (24, EMB_DIMS)
    assert np.isfinite(embeddings).all()
    assert output.mod["prot"].shape == (24, 2)

    # the same cells in reverse order must get the same embeddings, which only
    # holds when the length sort of EmbExtractor is undone
    reversed_h5mu = tokenized_h5mu.copy()
    reversed_h5mu.mod["rna"] = tokenized_h5mu.mod["rna"][::-1].copy()
    reversed_h5mu.mod["prot"] = tokenized_h5mu.mod["prot"][::-1].copy()
    reversed_path = random_h5mu_path()
    reversed_h5mu.write_h5mu(reversed_path)
    reversed_output = run_embed(run_component, reversed_path, random_h5mu_path())
    reversed_embeddings = reversed_output.mod["rna"].obsm["X_geneformer"]

    assert np.allclose(embeddings, reversed_embeddings[::-1], atol=1e-5)
    # and different cells get different embeddings
    assert not np.allclose(embeddings[0], embeddings[1])


def test_custom_obsm_output(run_component, input_path, random_h5mu_path):
    output = run_embed(
        run_component,
        input_path,
        random_h5mu_path(),
        "--obsm_output",
        "X_custom",
        "--emb_layer",
        "-1",
        "--output_compression",
        "gzip",
    )
    assert output.mod["rna"].obsm["X_custom"].shape == (24, EMB_DIMS)
    assert "X_geneformer" not in output.mod["rna"].obsm


def test_model_version_mismatch_errors(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--model",
                model_dir,
                "--model_version",
                "V2",
            ]
        )
    assert re.search(
        r"tokenized for model version V1, but --model_version is V2",
        err.value.stdout.decode("utf-8"),
    )


def test_cls_mode_needs_v2(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_embed(run_component, input_path, random_h5mu_path(), "--emb_mode", "cls")
    assert re.search(r"V1 models do not have", err.value.stdout.decode("utf-8"))


def test_untokenized_input_errors(run_component, tokenized_h5mu, random_h5mu_path):
    del tokenized_h5mu.mod["rna"].obsm["geneformer_tokens"]
    path = random_h5mu_path()
    tokenized_h5mu.write_h5mu(path)
    with pytest.raises(CalledProcessError) as err:
        run_embed(run_component, path, random_h5mu_path())
    assert re.search(
        r"'geneformer_tokens' not found in .obsm", err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
