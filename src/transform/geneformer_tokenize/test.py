import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from datasets import load_from_disk
from geneformer import TranscriptomeTokenizer

## VIASH START
meta = {
    "executable": "./target/executable/transform/geneformer_tokenize/geneformer_tokenize",
    "resources_dir": "resources_test",
    "config": "./src/transform/geneformer_tokenize/config.vsh.yaml",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def input_h5mu():
    # 60 cells keep the test fast; the gene space is left untouched because the
    # rank of a gene depends on all genes of the cell
    mdata = mu.read_h5mu(input_file)
    return mu.MuData({mod: mdata.mod[mod][:60].copy() for mod in mdata.mod})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write_h5mu(path)
    return path


def reference_tokens(adata, model_version, tmp_path):
    """Tokenize with TranscriptomeTokenizer.tokenize_data, the upstream entry
    point, and return the input_ids per cell in input order."""
    staged = ad.AnnData(
        X=adata.X,
        obs=pd.DataFrame(
            {
                "n_counts": np.asarray(adata.X.sum(axis=1)).ravel(),
                "row": np.arange(adata.n_obs),
            },
            index=adata.obs_names,
        ),
        var=pd.DataFrame(
            {"ensembl_id": adata.var_names.to_numpy()}, index=adata.var_names
        ),
    )
    in_dir = tmp_path / f"reference_in_{model_version}"
    in_dir.mkdir()
    staged.write_h5ad(in_dir / "cells.h5ad")
    tokenizer = TranscriptomeTokenizer(
        {"row": "row"}, nproc=1, model_version=model_version
    )
    tokenizer.tokenize_data(
        in_dir, tmp_path, f"reference_{model_version}", file_format="h5ad"
    )
    dataset = load_from_disk(str(tmp_path / f"reference_{model_version}.dataset"))
    by_row = dict(zip(dataset["row"], dataset["input_ids"]))
    return [np.asarray(by_row[i]) for i in range(adata.n_obs)], tokenizer


@pytest.mark.parametrize("model_version,width", [("V2", 4096), ("V1", 2048)])
def test_tokens_match_upstream_tokenizer(
    run_component,
    input_path,
    input_h5mu,
    random_h5mu_path,
    tmp_path,
    model_version,
    width,
):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--model_version",
            model_version,
            "--output_compression",
            "gzip",
        ]
    )

    output = mu.read_h5mu(output_path)
    rna = output.mod["rna"]
    tokens = rna.obsm["geneformer_tokens"]
    assert tokens.shape == (60, width)

    expected, tokenizer = reference_tokens(
        input_h5mu.mod["rna"], model_version, tmp_path
    )
    pad = tokenizer.gene_token_dict["<pad>"]
    for row, want in enumerate(expected):
        got = tokens[row]
        assert np.array_equal(got[: len(want)], want), (
            f"row {row} differs from upstream"
        )
        assert np.all(got[len(want) :] == pad), f"row {row} is not padded with <pad>"

    settings = rna.uns["geneformer_tokenize"]
    assert settings["model_version"] == model_version
    assert settings["model_input_size"] == width
    assert bool(settings["special_token"]) == (model_version == "V2")
    if model_version == "V2":
        assert np.all(tokens[:, 0] == settings["cls_token"])

    # every gene with a token is a gene the tokenizer knows
    known = rna.var["geneformer_token"] >= 0
    assert known.sum() > 10_000
    assert set(rna.var.loc[known, "geneformer_token"]) <= set(
        tokenizer.gene_token_dict.values()
    )
    assert (rna.var.loc[~known, "geneformer_ensembl_id"] == "").all()

    # other modalities are not touched
    assert output.mod["prot"].shape == input_h5mu.mod["prot"].shape


def test_uncropped_tokens_extend_model_input(
    run_component, input_path, random_h5mu_path
):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--model_version",
            "V1",
            "--obsm_output_uncropped",
            "geneformer_tokens_full",
        ]
    )

    rna = mu.read_h5mu(output_path).mod["rna"]
    cropped = rna.obsm["geneformer_tokens"]
    full = rna.obsm["geneformer_tokens_full"]
    pad = rna.uns["geneformer_tokenize"]["pad_token"]

    n_genes = (full != pad).sum(axis=1)
    n_expressed_known = np.array(
        [
            np.isin(
                rna.var["geneformer_token"].to_numpy()[row.indices], [-1], invert=True
            ).sum()
            for row in rna.X
        ]
    )
    assert np.array_equal(n_genes, n_expressed_known), (
        "uncropped tokens should hold every expressed gene in the dictionary"
    )
    # V1 has no special tokens, so the model input is the uncropped prefix
    width = min(cropped.shape[1], full.shape[1])
    assert np.array_equal(cropped[:, :width], full[:, :width])


def test_gene_map_gives_same_tokens(
    run_component, input_h5mu, random_h5mu_path, tmp_path
):
    # symbols as .var index, Ensembl ids only available through the gene map
    by_id_path = random_h5mu_path()
    input_h5mu.write_h5mu(by_id_path)

    rna = input_h5mu.mod["rna"].copy()
    ensembl_ids = rna.var_names.to_numpy()
    rna.var_names = rna.var["gene_symbol"].astype(str).to_numpy()
    rna.var_names_make_unique()
    gene_map = pd.DataFrame({"gene_symbol": rna.var_names, "ensembl_id": ensembl_ids})
    gene_map_path = tmp_path / "gene_map.csv"
    gene_map.to_csv(gene_map_path, index=False)
    by_symbol_path = random_h5mu_path()
    mu.MuData({"rna": rna}).write_h5mu(by_symbol_path)

    out_by_id = random_h5mu_path()
    out_by_symbol = random_h5mu_path()
    run_component(["--input", by_id_path, "--output", out_by_id])
    run_component(
        [
            "--input",
            by_symbol_path,
            "--output",
            out_by_symbol,
            "--gene_map",
            gene_map_path,
        ]
    )

    tokens_by_id = mu.read_h5mu(out_by_id).mod["rna"].obsm["geneformer_tokens"]
    tokens_by_symbol = mu.read_h5mu(out_by_symbol).mod["rna"].obsm["geneformer_tokens"]
    assert np.array_equal(tokens_by_id, tokens_by_symbol)


def test_normalized_counts_are_refused(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--input_layer",
                "log_normalized",
            ]
        )
    assert re.search(r"non-integer values", err.value.stdout.decode("utf-8"))

    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--input_layer",
            "log_normalized",
            "--allow_non_integer",
        ]
    )
    assert "geneformer_tokens" in mu.read_h5mu(output_path).mod["rna"].obsm


def test_missing_layer_errors(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--input_layer",
                "does_not_exist",
            ]
        )
    assert re.search(
        r"Layer 'does_not_exist' not found", err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
